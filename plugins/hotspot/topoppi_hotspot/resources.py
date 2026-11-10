"""Locate the bundled official ranking model or an explicitly supplied model."""
from importlib import resources
from pathlib import Path


def default_model_dir():
    """Return the bundled mutation-derived surface model."""
    return Path(str(resources.files("topoppi_hotspot").joinpath("data", "mutation_surface")))


def resolve_model_dir(model_dir=None):
    directory = Path(model_dir).expanduser() if model_dir else default_model_dir()
    if not (directory / "model.json").is_file():
        raise ValueError("Choose a model folder containing model.json, or leave it empty to use the official model.")
    return directory.resolve()
