"""Public model vocabulary preserves fitted state and scientific provenance."""
import json

import pytest

from topoppi_affinity.cli import main
from topoppi_affinity.model_metadata import training_metadata


def test_metadata_preserves_membership_and_numeric_protocol():
    raw = {
        "representation": "legacy_v1", "n_train": 2,
        "train_sample_ids": ["record001", "record002"],
        "parameters": {"alpha": 1.0, "legacy_weight": .25},
        "selection": "Fixed equal weights; no confirmation labels",
        "root": "/example/private-experiment", "confirmation_evaluated": False,
        "checkpoint": "/example/private-experiment/model.pt",
        "prior_provenance": {"method": "StaB-ddG-stability-only", "source_url": "https://example.org/model",
                             "training_data": "ProteinMPNN structures and Megascale folding stability",
                             "checkpoint": "/local/research/model_ckpts/stability_finetuned.pt", "mc_samples": 20},
    }
    result = training_metadata(raw)
    assert result["train_sample_ids"] == raw["train_sample_ids"]
    assert result["n_train"] == 2
    assert result["checkpoint"] == "model.pt"
    assert result["parameters"] == {"alpha": 1.0, "descriptor_kernel_weight": .25}
    assert result["representation"] == "conventional_descriptors"
    assert result["prior_provenance"]["checkpoint"] == "stability_finetuned.pt"
    assert result["prior_provenance"]["mc_samples"] == 20
    assert training_metadata(result) == result
    assert all(word not in json.dumps(result) for word in ("private-experiment", "confirmation", "legacy_v1"))








@pytest.mark.parametrize("arguments,expected", [
    (["predict", "--help"], "--descriptor-features"),
    (["review", "--help"], "--view {structure,contacts,atlas}"),
])
def test_public_cli_uses_task_and_view_names(arguments, expected, capsys):
    with pytest.raises(SystemExit) as exc:
        main(arguments)
    assert exc.value.code == 0
    text = capsys.readouterr().out
    assert expected in text
    assert "legacy" not in text.lower() and "--condition" not in text


def test_public_cli_does_not_dispatch_experiment_commands(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["run-pilots"])
    assert exc.value.code == 2
    assert "invalid choice" in capsys.readouterr().err
