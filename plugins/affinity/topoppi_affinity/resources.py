"""Explicit, resumable preparation of user-selected prediction resources.

Imports, status inspection and model resolution never access the network. ESMC
is fetched only from its official host, never from a TopoPPI distribution.
"""
from __future__ import annotations

import http.client
import json
import os
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from importlib.resources import files
from pathlib import Path

TASKS = ("absolute", "mutation", "known", "unseen")
RELEASE = "3.0.0"
ESMC_FILE = "esmc_600m_2024_12_v0.pth"


def _preferences_path():
    return Path.home() / ".topoppi" / "affinity.json"


def remember_cache(cache_dir):
    """Remember a user's explicit resource-location choice for later CLI/GUI runs."""
    path = Path(cache_dir).expanduser().resolve()
    config = _preferences_path()
    config.parent.mkdir(parents=True, exist_ok=True)
    config.write_text(json.dumps({"cache_dir": str(path)}, indent=2) + "\n", encoding="utf-8")
    return path


def cache_root(cache_dir=None):
    """Shared model/embedding cache; TOPOPPI_AFFINITY_CACHE overrides defaults."""
    if cache_dir is not None:
        return Path(cache_dir).expanduser()
    configured = os.environ.get("TOPOPPI_AFFINITY_CACHE")
    if configured:
        return Path(configured).expanduser()
    preferences = _preferences_path()
    if preferences.is_file():
        chosen = json.loads(preferences.read_text(encoding="utf-8")).get("cache_dir")
        if chosen:
            return Path(chosen).expanduser()
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData/Local"))
    elif sys.platform == "darwin":
        base = Path.home() / "Library/Caches"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base / "topoppi" / "affinity"


def model_catalog(catalog_path=None):
    source = Path(catalog_path) if catalog_path else files("topoppi_affinity").joinpath("data/model_catalog.json")
    return json.loads(source.read_text(encoding="utf-8"))


def resolve_model(task="absolute", model=None, *, cache_dir=None):
    if model is not None and str(model).strip():
        path = Path(model).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Model does not exist: {path}")
        return path.resolve()
    if task not in TASKS:
        raise ValueError(f"Select one of: {', '.join(TASKS)}")
    path = cache_root(cache_dir) / "models" / RELEASE / task / "model.joblib"
    if not path.is_file() and task in {"absolute", "mutation"}:
        from topoppi_affinity.production import bundled_model
        path = bundled_model(task)
    if not path.is_file():
        raise FileNotFoundError(
            f"Choose a local model for {task} with --model or Choose local model in the Affinity window. "
            f"To prepare a supplied catalog, run `topoppi-affinity models download --task {task} --catalog CATALOG.json`."
        )
    return path.resolve()


def esmc_checkpoint(*, cache_dir=None):
    """Locate our download or the same frozen weights already in the SDK cache."""
    prepared = cache_root(cache_dir) / "upstream" / "esmc" / ESMC_FILE
    if prepared.is_file():
        return prepared
    hf_home = Path(os.environ.get("HF_HOME", Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "huggingface"))
    hub = Path(os.environ.get("HF_HUB_CACHE", os.environ.get("HUGGINGFACE_HUB_CACHE", hf_home / "hub")))
    revision = model_catalog()["esmc"]["revision"]
    for owner in ("biohub", "EvolutionaryScale"):
        path = hub / f"models--{owner}--esmc-600m-2024-12" / "snapshots" / revision / "data/weights" / ESMC_FILE
        if path.is_file():
            return path
    return None


def prior_options(task, *, cache_dir=None):
    if task not in {"known", "unseen"}:
        return {}
    spec = model_catalog()["priors"][task]
    repository = cache_root(cache_dir) / "upstream" / spec["directory"]
    if not _prior_ready(repository, spec):
        return {}
    result = {"repository": repository}
    if spec.get("checkpoint"):
        result["checkpoint"] = repository / spec["checkpoint"]
    return result


def _prior_ready(repository, spec):
    # Downloaded BA-Cycle code has no .git directory: inference reads this
    # provenance record. It is written only after all resources are prepared.
    return ((repository / ".topoppi-source.json").is_file()
            and all((repository / name).is_file() for name in spec["required_files"]))


def resource_status(*, cache_dir=None):
    catalog = model_catalog()
    root = cache_root(cache_dir)
    models = {}
    checkpoint = esmc_checkpoint(cache_dir=cache_dir)
    for task, spec in catalog["models"].items():
        path = root / "models" / RELEASE / task / "model.joblib"
        bundled = False
        if not path.is_file() and task in {"absolute", "mutation"}:
            from topoppi_affinity.production import bundled_model
            path = bundled_model(task)
            bundled = path.is_file()
        prior = prior_options(task, cache_dir=cache_dir)
        missing = []
        if not path.is_file():
            missing.append(spec["name"] + " model")
        if not bundled and task in {"known", "unseen"} and not prior:
            missing.append(catalog["priors"][task]["name"])
        if not bundled and checkpoint is None:
            missing.append("ESMC-600M sequence weights")
        models[task] = {"ready": not missing, "missing_resources": missing, "name": spec["name"], "target": spec["target"], "scope": spec["scope"],
            "path": str(path), "installed": path.is_file(), "bundled": bundled, "prior_required": task in {"known", "unseen"},
            "prior_ready": task not in {"known", "unseen"} or bool(prior),
            "prior_license": catalog["priors"].get(task, {}).get("license"),
            "prior": {key: str(value) for key, value in prior.items()}, "url": spec["url"]}
    return {"release": RELEASE, "cache_dir": str(root), "models": models,
        "esmc": {"installed": checkpoint is not None, "path": str(checkpoint) if checkpoint else None,
                 "url": catalog["esmc"]["url"], "size_bytes": catalog["esmc"]["size_bytes"],
                 "source": "Downloaded directly from the official ESMC repository"}}


def _check_cancel(cancel_event):
    if cancel_event is not None and cancel_event.is_set():
        raise InterruptedError("Resource preparation cancelled; downloaded partial files are kept for retry")


def _notify(progress, event, resource, current=0, total=None, message=None):
    if progress:
        progress({"event": event, "resource": resource, "current": current, "total": total,
                  "message": message or f"{resource}: {current / 1_000_000:.1f} MB"
                  + (f" / {total / 1_000_000:.1f} MB" if total else "")})


def _download(url, destination, *, resource, progress=None, cancel_event=None, expected_size=None):
    """Keep interrupted transfers and resume when the official host supports Range."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file():
        if expected_size is None or destination.stat().st_size == expected_size:
            _notify(progress, "cached", resource, message=f"{resource}: already downloaded")
            return destination
        destination.unlink()
    partial = destination.with_name(destination.name + ".part")
    if expected_size is not None:
        remaining = max(0, expected_size - (partial.stat().st_size if partial.exists() else 0))
        if shutil.disk_usage(destination.parent).free < remaining:
            raise OSError(f"Not enough disk space for {resource}: {remaining / 1e9:.2f} GB more is required. "
                          "Free space or choose another --cache-dir.")
    for attempt in range(3):
        _check_cancel(cancel_event)
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": "TopoPPI/3.0 affinity-resource-preparation"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=30) as response:
                resumed = offset > 0 and response.status == 206
                if not resumed:
                    offset = 0
                length = response.headers.get("Content-Length")
                total = (offset + int(length)) if length else expected_size
                _notify(progress, "downloading", resource, offset, total)
                with partial.open("ab" if resumed else "wb") as handle:
                    current, last_update = offset, 0.
                    while True:
                        _check_cancel(cancel_event)
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                        current += len(chunk)
                        if time.monotonic() - last_update >= .2:
                            _notify(progress, "downloading", resource, current, total)
                            last_update = time.monotonic()
                if (total is not None and current != total) or (expected_size is not None and current != expected_size):
                    raise OSError(f"Incomplete download of {resource}: {current} bytes; expected {expected_size or total}")
            partial.replace(destination)
            _notify(progress, "downloaded", resource, current, total)
            return destination
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and partial.exists():
                if expected_size and partial.stat().st_size == expected_size:
                    partial.replace(destination)
                    return destination
                partial.unlink()
            elif exc.code == 404:
                raise OSError(f"{resource} is not published at {url}. Use the matching published TopoPPI release "
                              "or provide your own compatible model with --model.") from exc
            elif exc.code not in {408, 429, 500, 502, 503, 504}:
                raise OSError(f"Could not download {resource}: HTTP {exc.code} ({url})") from exc
            if attempt == 2:
                raise OSError(f"Could not finish {resource}; run preparation again to retry: {exc}") from exc
        except InterruptedError:
            raise
        except (OSError, urllib.error.URLError, http.client.HTTPException) as exc:
            if attempt == 2:
                raise OSError(f"Could not finish {resource}; run preparation again to resume: {exc}") from exc
        _notify(progress, "retry", resource, message=f"{resource}: retrying interrupted download ({attempt + 2}/3)")
    return destination


def _extract(archive, destination, *, strip_root=False, cancel_event=None):
    """Extract published archives into a sibling staging directory before use."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".preparing-", dir=destination.parent) as temporary:
        staging = Path(temporary)
        with zipfile.ZipFile(archive) as source:
            for member in source.infolist():
                _check_cancel(cancel_event)
                parts = Path(member.filename).parts[1:] if strip_root else Path(member.filename).parts
                if not parts or member.is_dir():
                    continue
                if ".." in parts or Path(*parts).is_absolute():
                    raise ValueError(f"Invalid archive path: {member.filename}")
                target = staging.joinpath(*parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                with source.open(member) as read, target.open("wb") as write:
                    shutil.copyfileobj(read, write)
        if destination.exists():
            shutil.rmtree(destination)
        staging.rename(destination)


def _download_prior(task, catalog, root, *, progress, cancel_event):
    spec = catalog["priors"][task]
    destination = root / "upstream" / spec["directory"]
    if _prior_ready(destination, spec):
        _notify(progress, "cached", spec["name"], message=f"{spec['name']}: ready")
        return
    _notify(progress, "source", spec["name"], message=f"{spec['name']}: {spec['source_url']} ({spec['license']})")
    archive = _download(spec["url"], root / "downloads" / (spec["directory"] + ".zip"),
        resource=spec["name"], progress=progress, cancel_event=cancel_event)
    _extract(archive, destination, strip_root=True, cancel_event=cancel_event)
    for item in spec.get("downloads", []):
        prepared = _download(item["url"], root / "downloads" / (spec["directory"] + "-" + Path(item["path"]).name), resource=item["name"],
            progress=progress, cancel_event=cancel_event, expected_size=item.get("size_bytes"))
        target = destination / item["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(prepared, target)
    if not all((destination / name).is_file() for name in spec["required_files"]):
        raise OSError(f"{spec['name']} preparation is incomplete; rerun preparation to retry")
    (destination / ".topoppi-source.json").write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    _notify(progress, "ready", spec["name"], message=f"{spec['name']}: ready")


def download_models(tasks=("absolute",), *, include_esmc=True, include_priors=True, cache_dir=None, progress=None, cancel_event=None,
                    catalog_path=None):
    """Prepare selected frozen models and their priors; called only by user action."""
    tasks = tuple(dict.fromkeys((tasks,) if isinstance(tasks, str) else tasks))
    if (not tasks and not include_esmc) or any(task not in TASKS for task in tasks):
        raise ValueError(f"Select at least one task from {', '.join(TASKS)}")
    _check_cancel(cancel_event)
    if cache_dir is not None:
        remember_cache(cache_dir)
    catalog, root = model_catalog(catalog_path), cache_root(cache_dir)
    bundled_only = bool(tasks) and catalog_path is None and all(t in {"absolute", "mutation"} for t in tasks)
    for task in tasks:
        _check_cancel(cancel_event)
        if catalog_path is None and task in {"absolute", "mutation"}:
            _notify(progress, "ready", task, message=f"{task}: bundled default model ready")
            continue
        spec = catalog["models"][task]
        destination = root / "models" / RELEASE / task
        if not (destination / "model.joblib").is_file():
            if not spec.get("url"):
                raise ValueError(
                    f"No model source is configured for {task}. Supply --catalog with an explicit "
                    "model archive URL, or use --model with a compatible local fitted model."
                )
            archive = _download(spec["url"], root / "downloads" / Path(spec["url"]).name,
                resource=spec["name"], progress=progress, cancel_event=cancel_event,
                expected_size=spec.get("size_bytes"))
            _extract(archive, destination, cancel_event=cancel_event)
            if not (destination / "model.joblib").is_file():
                raise ValueError(f"The {task} release archive contains no model.joblib")
        _notify(progress, "ready", task, message=f"{spec['name']}: ready")
        if include_priors and task in catalog["priors"]:
            _download_prior(task, catalog, root, progress=progress, cancel_event=cancel_event)
    if include_esmc and not bundled_only:
        checkpoint = esmc_checkpoint(cache_dir=cache_dir)
        if checkpoint:
            _notify(progress, "cached", "ESMC-600M", message=f"ESMC-600M: reusing {checkpoint}")
        else:
            spec = catalog["esmc"]
            _download(spec["url"], root / "upstream/esmc" / ESMC_FILE, resource="ESMC-600M (official upstream)",
                progress=progress, cancel_event=cancel_event, expected_size=spec.get("size_bytes"))
    return resource_status(cache_dir=cache_dir)


def export_example(output_dir, *, reuse=False):
    """Copy the installed, successful 1BVK example without needing a git checkout."""
    destination = Path(output_dir)
    source = files("topoppi_affinity").joinpath("data/example")
    names = ("1bvk.pdb", "variants.csv", "README.md")
    existing = [name for name in names if (destination / name).exists()]
    if existing and (not reuse or any((destination / name).read_bytes() != source.joinpath(name).read_bytes() for name in existing)):
        raise FileExistsError(f"Example files already exist or have been edited in {destination}; choose a new directory or analyze those files explicitly")
    destination.mkdir(parents=True, exist_ok=True)
    for name in names:
        (destination / name).write_bytes(source.joinpath(name).read_bytes())
    return {"structure": str(destination / "1bvk.pdb"), "variants": str(destination / "variants.csv"),
            "chains_a": "A,B", "chains_b": "C", "task": "absolute"}
