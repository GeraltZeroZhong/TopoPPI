"""Small file operations shared across TopoPPI entry points."""

from __future__ import annotations

import csv
import hashlib
import os
import subprocess
import tempfile
from os import PathLike
from pathlib import Path
from typing import Mapping, Sequence


def sha256_file(path: str | PathLike[str]) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_worktree_state(repo_root: str | PathLike[str]) -> tuple[str | None, bool | None]:
    """Return revision metadata when *repo_root* is a Git worktree."""

    root = Path(repo_root)
    if not (root / ".git").exists():
        return None, None
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None, None
    revision = commit.stdout.strip() if commit.returncode == 0 else None
    dirty = bool(status.stdout.strip()) if status.returncode == 0 else None
    return revision, dirty


def read_csv_rows(path: str | PathLike[str]) -> list[dict[str, str]]:
    """Read a header-based UTF-8 CSV file."""

    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv_atomic(
    path: str | PathLike[str],
    rows: list[Mapping[str, object]],
    fieldnames: Sequence[str] | None = None,
) -> None:
    """Write CSV rows beside their destination, then replace it."""

    target = Path(path)
    fields = list(fieldnames) if fieldnames is not None else sorted({key for row in rows for key in row})
    with tempfile.NamedTemporaryFile("w", newline="", encoding="utf-8", dir=target.parent, delete=False) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, target)


def write_npz_atomic(path: str | PathLike[str], arrays: Mapping[str, object]) -> None:
    """Replace a compressed NPZ only after the complete archive has been written."""
    import numpy as np

    target = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("wb", dir=target.parent, prefix=f".{target.name}.",
                                         suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            np.savez_compressed(handle, **arrays)
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


__all__ = ["git_worktree_state", "read_csv_rows", "sha256_file", "write_csv_atomic", "write_npz_atomic"]
