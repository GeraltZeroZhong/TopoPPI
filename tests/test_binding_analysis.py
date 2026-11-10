"""Scientific identities and explicit grouping survive cross-task navigation."""
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from topoppi.gui_app.binding_analysis import (
    BindingAnalysisMixin,
    normalize_residues,
    validate_variant_structure,
    variant_rows,
)


def test_variant_grouping_preserves_chains_author_numbers_and_insertion_codes():
    residues = normalize_residues([
        {"chain": "B", "position": 52, "insertion": "", "wt": "Y"},
        {"chain": "B", "position": 52, "insertion": "A", "wt": "W"},
        {"chain": "C", "position": 52, "insertion": "", "wt": "G"},
        {"chain": "B", "position": 52, "insertion": "", "wt": "Y"},
    ])
    targets = dict.fromkeys(residues, "A")
    singles = variant_rows(residues, targets, "single")
    assert [row["mutations"] for row in singles] == ["B:Y52A", "B:W52AA", "C:G52A"]
    assert variant_rows(residues, targets, "joint")[0]["mutations"] == "B:Y52A,B:W52AA,C:G52A"


def test_native_alanine_needs_an_explicit_non_identity_target():
    residues = normalize_residues([{"chain": "A", "position": -1, "insertion": "", "wt": "A"}])
    with pytest.raises(ValueError, match="already A"):
        variant_rows(residues, {"A:-1": "A"}, "single")
    assert variant_rows(residues, {"A:-1": "V"}, "single")[0]["mutations"] == "A:A-1V"


def test_unknown_wildtype_and_conflicting_identity_are_not_guessed():
    with pytest.raises(ValueError, match="standard amino acid"):
        normalize_residues([{"chain": "A", "position": 3, "wt": "X"}])
    with pytest.raises(ValueError, match="Conflicting identities"):
        normalize_residues([{"chain": "A", "position": 3, "wt": "G"}, {"chain": "A", "position": 3, "wt": "Y"}])


def test_navigation_passes_explicit_task_and_source_without_running_predictions():
    host = BindingAnalysisMixin()
    callback = Mock(return_value="panel")
    host._binding_tasks = {}
    host._busy = host._closing = False
    host.register_binding_task("mutation", title="Mutation effects", description="", plugin="affinity", open_task=callback)
    setup = {"structure": "original.pdb", "chains_a": "A,B", "chains_b": "C", "append_variants": [{"variant_id": "x", "mutations": "B:G42V"}]}
    assert host.open_binding_task("mutation", setup) == "panel"
    assert callback.call_args.kwargs["setup"] == {**setup, "prediction_task": "mutation"}
    assert "prediction_task" not in setup
    host._busy = True
    with pytest.raises(RuntimeError, match="Wait"):
        host.open_binding_task("mutation", setup)
    assert callback.call_count == 1


def test_missing_plugin_and_missing_coordinate_source_have_actionable_recovery(tmp_path):
    host = BindingAnalysisMixin()
    host._binding_tasks = {}
    host._busy = host._closing = False
    with pytest.raises(RuntimeError, match="Open Plugins"):
        host.open_binding_task("mutation", {})
    atlas = tmp_path / "saved.npz"
    atlas.touch()
    with pytest.raises(ValueError, match="existing PDB or mmCIF"):
        validate_variant_structure(atlas, {}, {})


def test_relocated_variant_source_must_match_recorded_partners_and_residue_identity():
    from pathlib import Path

    path = Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"
    rows = normalize_residues([{"chain": "A", "position": 1, "wt": "A"}])
    assert validate_variant_structure(path, {"chains_a": "A", "chains_b": "B"}, rows) == str(path.resolve())
    with pytest.raises(ValueError, match="partner chains"):
        validate_variant_structure(path, {"chains_a": "A", "chains_b": "C"}, rows)
    rows["A:1"]["wt"] = "G"
    with pytest.raises(ValueError, match="does not match selected residue A:G1"):
        validate_variant_structure(path, {"chains_a": "A", "chains_b": "B"}, rows)


def test_unreadable_relocated_structure_has_recovery_message(tmp_path):
    path = tmp_path / "incomplete.cif"
    path.write_text("data_incomplete\n")
    with pytest.raises(ValueError, match="Choose a readable PDB/mmCIF"):
        validate_variant_structure(path, {"chains_a": "A", "chains_b": "B"}, {})


def test_displayed_source_does_not_follow_edited_input_fields():
    host = BindingAnalysisMixin()
    host.var_input_path = SimpleNamespace(get=lambda: "next.pdb")
    host.var_chain_a = SimpleNamespace(get=lambda: "X")
    host.var_chain_b = SimpleNamespace(get=lambda: "Y")
    host._successful_single_run = {"params": {"path": "original.pdb", "chain_a": "A,B", "chain_b": "C"},
                                   "viz": SimpleNamespace(chain_a_id="A,B", chain_b_id="C")}
    assert host.binding_setup(displayed=True) == {"structure": "original.pdb", "chains_a": "A,B", "chains_b": "C"}
    assert host.binding_setup()["structure"] == "next.pdb"


def test_core_selection_uses_the_same_modified_residue_identity_as_prediction():
    host = BindingAnalysisMixin()
    host._busy = False
    host._binding_selection = {}
    host._show_binding_selection = Mock()
    host._successful_single_run = {"viz": SimpleNamespace(residue_metadata_B={}, residue_metadata_A={
        "A:MSE:5": {"residue_name": "MSE", "residue_id": ("H_MSE", 5, " "), "chain": "A"}})}
    host._toggle_binding_residue("A:MSE:5")
    assert host._binding_selection["A:MSE:5"] == {"chain": "A", "position": 5, "insertion": "", "wt": "M"}
    host._toggle_binding_residue("A:MSE:5")
    assert not host._binding_selection


def test_selection_marker_uses_the_residues_own_patch_axis():
    host = BindingAnalysisMixin()
    host._busy = False
    host._binding_selection = {"A:TYR:2": {}}
    host._binding_selection_artists = []
    host.binding_selection_label = Mock()
    first, second = Mock(), Mock()
    host.current_fig = SimpleNamespace(axes=[first, second])
    host.current_canvas = Mock()
    host._surface_view = lambda: None
    host._label_anchor = lambda objects: objects["anchor"]
    host._successful_single_run = {"viz": SimpleNamespace(artist_map={"patch2": {
        "residue_key": "A:TYR:2", "anchor": [2.0, 3.0], "collection": SimpleNamespace(axes=second)}})}
    host._show_binding_selection()
    first.scatter.assert_not_called()
    second.scatter.assert_called_once()
