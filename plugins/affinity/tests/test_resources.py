"""Exercise the installation/retry paths used by CLI and GUI, without Internet."""
from __future__ import annotations

import http.client
import io
import json
import threading
import zipfile

import pytest

from topoppi_affinity import resources


@pytest.fixture(autouse=True)
def isolated_preferences(tmp_path, monkeypatch):
    monkeypatch.setattr(resources, "_preferences_path", lambda: tmp_path / "preferences.json")


def test_status_and_resolution_are_offline(tmp_path, monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("A read-only operation attempted a download")

    monkeypatch.setattr(resources.urllib.request, "urlopen", unexpected)
    monkeypatch.setenv("TOPOPPI_AFFINITY_CACHE", str(tmp_path))
    assert resources.cache_root() == tmp_path
    assert resources.resource_status()["models"]["absolute"]["installed"]
    assert resources.resolve_model().name == "model.json"
    custom = tmp_path / "custom.joblib"
    custom.write_bytes(b"existing custom model")
    assert resources.resolve_model("absolute", custom) == custom


def test_local_release_install_and_repeat_skips_network(tmp_path, monkeypatch):
    archive = tmp_path / "absolute.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("model.joblib", b"frozen model bytes")
        handle.writestr("model_report.json", '{"target":"pkd"}')
    catalog = resources.model_catalog()
    catalog["models"]["absolute"].update(url=archive.as_uri(), size_bytes=archive.stat().st_size)
    catalog_path = tmp_path / "catalog.json"
    catalog_path.write_text(json.dumps(catalog))
    cache = tmp_path / "cache"
    updates = []
    state = resources.download_models("absolute", include_esmc=False, cache_dir=cache,
        catalog_path=catalog_path, progress=updates.append)
    assert state["models"]["absolute"]["installed"]
    assert resources.resolve_model(cache_dir=cache).read_bytes() == b"frozen model bytes"
    assert updates[-1]["event"] == "ready"
    assert resources.cache_root() == cache
    monkeypatch.setattr(resources.urllib.request, "urlopen", lambda *a, **kw: pytest.fail("Repeated download"))
    resources.download_models("absolute", include_esmc=False, cache_dir=cache, catalog_path=catalog_path)


def test_cancelled_download_resumes_from_partial_bytes(tmp_path, monkeypatch):
    payload = b"synthetic download bytes" * 250_000
    ranges = []

    class Response(io.BytesIO):
        def __init__(self, data, *, start=0):
            super().__init__(data[start:])
            self.status = 206 if start else 200
            self.headers = {"Content-Length": str(len(data) - start)}

    def urlopen(request, **kwargs):
        requested = request.get_header("Range")
        ranges.append(requested)
        start = int(requested.split("=")[1].split("-")[0]) if requested else 0
        return Response(payload, start=start)

    monkeypatch.setattr(resources.urllib.request, "urlopen", urlopen)
    target = tmp_path / "weights.pth"
    cancel = threading.Event()

    def stop_after_first_chunk(record):
        if record["current"] > 0:
            cancel.set()

    with pytest.raises(InterruptedError, match="cancelled"):
        resources._download("https://example.invalid/weights", target, resource="weights",
            progress=stop_after_first_chunk, cancel_event=cancel, expected_size=len(payload))
    assert not target.exists()
    partial_size = target.with_name("weights.pth.part").stat().st_size
    assert 0 < partial_size < len(payload)
    cancel.clear()
    resources._download("https://example.invalid/weights", target, resource="weights",
        cancel_event=cancel, expected_size=len(payload))
    assert ranges[-1] == f"bytes={partial_size}-"
    assert target.read_bytes() == payload


def test_prior_options_need_all_inputs_and_example_is_self_contained(tmp_path):
    catalog = resources.model_catalog()
    assert resources.prior_options("unseen", cache_dir=tmp_path) == {}
    spec = catalog["priors"]["unseen"]
    repo = tmp_path / "upstream" / spec["directory"]
    for name in spec["required_files"]:
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("prepared")
    assert resources.prior_options("unseen", cache_dir=tmp_path) == {}
    (repo / '.topoppi-source.json').write_text(json.dumps(spec))
    assert resources.prior_options("unseen", cache_dir=tmp_path) == {
        "repository": repo, "checkpoint": repo / spec["checkpoint"]}
    exported = resources.export_example(tmp_path / "example")
    assert exported["chains_a"] == "A,B"
    assert "YA50F" in (tmp_path / "example/variants.csv").read_text()
    assert "ATOM" in (tmp_path / "example/1bvk.pdb").read_text()
    with pytest.raises(FileExistsError):
        resources.export_example(tmp_path / "example")


def test_esmc_and_priors_are_official_pinned_sources():
    catalog = resources.model_catalog()
    assert catalog["esmc"]["url"].startswith("https://huggingface.co/biohub/")
    assert catalog["esmc"]["revision"] in catalog["esmc"]["url"]
    for prior in catalog["priors"].values():
        assert prior["url"].startswith("https://codeload.github.com/")
        assert prior["commit"] in prior["url"]


def test_cache_preference_precedence_and_status_does_not_change_it(tmp_path, monkeypatch):
    chosen, env, explicit = (tmp_path / name for name in ("chosen", "env", "explicit"))
    resources.remember_cache(chosen)
    assert resources.cache_root() == chosen
    resources.resource_status(cache_dir=explicit)
    assert resources.cache_root() == chosen
    monkeypatch.setenv("TOPOPPI_AFFINITY_CACHE", str(env))
    assert resources.cache_root() == env
    assert resources.cache_root(explicit) == explicit


def test_esmc_only_preparation_reuses_official_cache(tmp_path, monkeypatch):
    prepared = tmp_path / 'official.pth'
    prepared.write_bytes(b'cached frozen checkpoint')
    monkeypatch.setattr(resources, 'esmc_checkpoint', lambda **kwargs: prepared)
    monkeypatch.setattr(resources.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Unexpected network access'))
    status = resources.download_models(tasks=(), include_esmc=True)
    assert status['esmc']['installed']
    assert all(status['models'][task]['ready'] for task in ('absolute', 'mutation'))
    assert not any(status['models'][task]['installed'] for task in ('known', 'unseen'))


def test_chunked_transfer_error_is_retried(tmp_path, monkeypatch):
    attempts = []
    payload = b'complete source archive'

    class Response(io.BytesIO):
        status = 200
        headers = {}

        def read(self, size=-1):
            if len(attempts) == 1:
                raise http.client.IncompleteRead(b'interrupted chunk')
            return super().read(size)

    def urlopen(*args, **kwargs):
        attempts.append(1)
        return Response(payload)

    monkeypatch.setattr(resources.urllib.request, 'urlopen', urlopen)
    path = resources._download('https://official.test/source.zip', tmp_path / 'source.zip', resource='source')
    assert len(attempts) == 2
    assert path.read_bytes() == payload


def test_cancelled_prior_metadata_download_is_completed_on_retry(tmp_path):
    catalog = resources.model_catalog()
    spec = catalog['priors']['unseen']
    source = tmp_path / 'source.zip'
    with zipfile.ZipFile(source, 'w') as archive:
        for name in spec['required_files']:
            if name != spec['checkpoint']:
                archive.writestr('upstream/' + name, b'upstream source')
    checkpoint = tmp_path / 'checkpoint.pt'
    checkpoint.write_bytes(b'frozen checkpoint')
    license_file = tmp_path / 'LICENSE'
    license_file.write_text('Official license')
    spec['url'] = source.as_uri()
    spec['downloads'][0].update(url=checkpoint.as_uri(), size_bytes=checkpoint.stat().st_size)
    spec['downloads'][1]['url'] = license_file.as_uri()
    cache, cancel = tmp_path / 'cache', threading.Event()

    def cancel_before_license(record):
        if record['event'] == 'downloaded' and record['resource'] == spec['downloads'][0]['name']:
            cancel.set()

    with pytest.raises(InterruptedError):
        resources._download_prior('unseen', catalog, cache, progress=cancel_before_license, cancel_event=cancel)
    assert (cache / 'upstream' / spec['directory'] / spec['checkpoint']).is_file()
    assert resources.prior_options('unseen', cache_dir=cache) == {}
    cancel.clear()
    resources._download_prior('unseen', catalog, cache, progress=None, cancel_event=cancel)
    assert resources.prior_options('unseen', cache_dir=cache)['checkpoint'].read_bytes() == b'frozen checkpoint'
    assert (cache / 'upstream' / spec['directory'] / '.topoppi-source.json').is_file()


def test_task_readiness_includes_weights_and_required_prior(tmp_path, monkeypatch):
    monkeypatch.setattr(resources, 'esmc_checkpoint', lambda **kwargs: None)
    cache = tmp_path / 'cache'
    model = cache / 'models' / resources.RELEASE / 'absolute' / 'model.joblib'
    model.parent.mkdir(parents=True)
    model.write_bytes(b'prepared model')
    state = resources.resource_status(cache_dir=cache)
    assert state['models']['absolute']['installed']
    assert not state['models']['absolute']['ready']
    assert state['models']['absolute']['missing_resources'] == ['ESMC-600M sequence weights']
    monkeypatch.setattr(resources, 'esmc_checkpoint', lambda **kwargs: tmp_path / 'weights.pth')
    state = resources.resource_status(cache_dir=cache)
    assert state['models']['absolute']['ready']
    assert not state['models']['known']['ready']
    assert resources.model_catalog()['priors']['known']['name'] in state['models']['known']['missing_resources']


@pytest.mark.parametrize('options,expected', [([], ('absolute',)),
    (['--all'], ('absolute', 'mutation')),
    (['--task', 'known', '--task', 'unseen'], ('known', 'unseen'))])
def test_cli_download_prepares_selected_tasks_only(tmp_path, monkeypatch, capsys, options, expected):
    from topoppi_affinity import cli
    calls = []
    state = resources.resource_status(cache_dir=tmp_path)
    def prepare(**kwargs):
        calls.append(kwargs)
        return state
    monkeypatch.setattr(resources, 'download_models', prepare)
    assert cli.main(['models', 'download', *options, '--cache-dir', str(tmp_path)]) == 0
    assert calls[0]['tasks'] == expected
    output = capsys.readouterr().out
    assert 'Model cache:' in output and ('Prepare:' in output or ': ready' in output)


def test_cli_status_offers_readable_and_json_output(tmp_path, monkeypatch, capsys):
    from topoppi_affinity import cli
    monkeypatch.setattr(resources.urllib.request, 'urlopen', lambda *a, **kw: pytest.fail('Status downloaded a file'))
    args = ['models', 'status', '--cache-dir', str(tmp_path)]
    assert cli.main(args) == 0
    output = capsys.readouterr().out
    assert str(tmp_path) in output and 'Absolute binding affinity: ready' in output
    assert cli.main([*args, '--json']) == 0
    state = json.loads(capsys.readouterr().out)
    assert state['cache_dir'] == str(tmp_path)
    assert state['models']['absolute']['ready']


def test_custom_prior_preparation_skips_official_prior(tmp_path, monkeypatch):
    cache = tmp_path / 'cache'
    model = cache / 'models' / resources.RELEASE / 'known' / 'model.joblib'
    model.parent.mkdir(parents=True)
    model.write_bytes(b'prepared model')
    monkeypatch.setattr(resources, '_download_prior', lambda *a, **kw: pytest.fail('Official prior requested despite override'))
    resources.download_models('known', include_esmc=False, include_priors=False, cache_dir=cache)


def test_unconfigured_task_source_requires_explicit_catalog(tmp_path, monkeypatch):
    monkeypatch.setattr(resources.urllib.request, "urlopen", lambda *a, **kw: pytest.fail("Unexpected network access"))
    with pytest.raises(ValueError, match="Supply --catalog"):
        resources.download_models("known", cache_dir=tmp_path, include_esmc=False)
