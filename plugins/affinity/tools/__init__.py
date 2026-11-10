"""Installed feature extraction and explicit model-fitting commands."""
from pathlib import Path


def tool_file(name):
    return str(Path(__file__).parent / name)
