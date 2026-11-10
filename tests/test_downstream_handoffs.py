"""Synthetic desktop journeys across prediction, ranking and editable atlases."""

import time
import tkinter as tk
from pathlib import Path

import pytest

pytest.importorskip("topoppi_affinity")
pytest.importorskip("topoppi_hotspot")

from test_atlas_io import make_atlas
from topoppi_affinity.workflow import analyze
from topoppi_hotspot.gui_workflow import RankingSettings, run_ranking

import topoppi
from topoppi.gui_app.app import InterfaceApp
from topoppi.visualization.atlas_io import save_atlas
from topoppi.visualization.visualizer import InterfaceVisualizer


def _button(parent, text):
    for child in parent.winfo_children():
        if child.winfo_class() in {"TButton", "Button"} and child.cget("text") == text:
            return child
        found = _button(child, text)
        if found is not None:
            return found
    return None


def _wait(root, panel):
    deadline = time.monotonic() + 60
    while panel._working and time.monotonic() < deadline:
        root.update()
        time.sleep(.02)
    root.update()
    assert not panel._working


@pytest.fixture
def desktop(tmp_path, monkeypatch):
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("A graphical display is required")
    root.withdraw()
    monkeypatch.setattr(InterfaceApp, "_recent_store_path", lambda _: tmp_path / "recent.json")
    errors = []
    monkeypatch.setattr("tkinter.messagebox.showerror", lambda *args, **kwargs: errors.append(args))
    app = InterfaceApp(root)
    root.update()
    try:
        if not {"absolute", "mutation", "hotspot"}.issubset(app._binding_tasks):
            pytest.skip("Both optional prediction plugins must be installed")
        yield app, errors
    finally:
        app.close()
        try:
            root.destroy()
        except tk.TclError:
            pass


def _synthetic_sources(tmp_path):
    """Extend the independently generated peptide fixture with author identities."""
    fixture = (Path(__file__).parent / "fixtures" / "synthetic_complex.pdb").read_text()
    source = tmp_path / "inserted.pdb"
    source.write_text("".join(
        line[:26] + "A" + line[27:]
        if line.startswith("ATOM  ") and line[21] == "A" and int(line[22:26]) == 2 else line
        for line in fixture.splitlines(keepends=True)
    ))
    # Three translated copies of the artificial chain form an independent
    # multichain input. No experimental sequence or coordinates are used.
    multichain = tmp_path / "three_chains.pdb"
    lines = ["REMARK Synthetic translated peptide assembly.\n"]
    serial = 1
    for chain, z_offset in (("A", 0.), ("B", 4.), ("C", 8.)):
        for line in fixture.splitlines():
            if line.startswith("ATOM  ") and line[21] == "A":
                z = float(line[46:54]) + z_offset
                lines.append(f"{line[:6]}{serial:5d}{line[11:21]}{chain}{line[22:46]}{z:8.3f}{line[54:]}\n")
                serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    multichain.write_text("".join(lines))
    return source, multichain


def test_navigation_uses_edited_inputs_but_ranked_variants_keep_result_source(desktop, tmp_path):
    app, errors = desktop
    source, multichain = _synthetic_sources(tmp_path)
    setup = {"structure": str(source), "chains_a": "A", "chains_b": "B"}
    analyze(str(source), "A", "B", output=str(tmp_path / "affinity"), task="absolute", device="cpu", log=lambda _: None)
    affinity = app.open_binding_task("absolute", setup)
    affinity.window.withdraw()
    affinity.load(tmp_path / "affinity" / "analysis.json")
    ranking = run_ranking(RankingSettings(str(source), "A", "B", "A", "", str(tmp_path / "ranking")))
    hotspot = app.open_binding_task("hotspot", setup)
    hotspot.window.withdraw()
    hotspot.display(ranking)

    next_setup = {"structure": str(multichain), "chains_a": "A,B", "chains_b": "C"}
    for name, value in next_setup.items():
        affinity.vars[name].set(value)
    _button(affinity.window, "Rank hotspots").invoke()
    app.root.update()
    assert hotspot.binding_setup() == next_setup
    _button(hotspot.window, "Binding affinity").invoke()
    app.root.update()
    assert affinity.binding_setup() == next_setup
    assert affinity.document.report["input"]["path"] == str(source)
    assert hotspot.document.metadata["structure"] == str(source)

    hotspot.table.selection_set("A:PHE:2A")
    hotspot.table.event_generate("<<TreeviewSelect>>")
    app.root.update()
    composer = hotspot.create_variants()
    assert composer.setup == setup
    assert composer.preview.item(composer.preview.get_children()[0], "values")[1] == "A:F2AA"
    composer.window.destroy()
    assert not errors


def test_main_selection_keeps_inserted_identity_after_editing_next_input_and_reopening(desktop, tmp_path):
    app, errors = desktop
    source, multichain = _synthetic_sources(tmp_path)
    setup = {"structure": str(source), "chains_a": "A", "chains_b": "B"}
    ranking = run_ranking(RankingSettings(str(source), "A", "B", "A", "", str(tmp_path / "ranking")))
    hotspot = app.open_binding_task("hotspot", setup)
    hotspot.window.withdraw()
    hotspot.display(ranking)
    hotspot.table.selection_set("A:PHE:2A")
    assert hotspot.show_surface()
    app.var_input_path.set(str(multichain))
    app.var_chain_a.set("A,B")
    app.var_chain_b.set("C")
    app._toggle_binding_residue("A:PHE:2A")
    composer = app._create_selected_variants()
    assert composer.setup == setup
    composer.add_button.invoke()
    affinity = app._affinity_panel
    affinity.window.withdraw()
    assert affinity.binding_setup() == setup
    assert affinity.variants.get("1.0", "end-1c") == "A_F2AA\tA:F2AA"
    affinity.vars["output"].set(str(tmp_path / "mutation"))
    affinity.run_button.invoke()
    _wait(app.root, affinity)
    assert affinity.document is not None, affinity.errors.get()
    assert affinity.document.predictions.status.eq("ok").all()
    assert app.affinity_table.heading("value", "text") == "ddG (kcal/mol)"

    affinity.load(affinity.document.path)
    affinity.show_in_main()
    app.root.update()
    app._clear_binding_selection()
    app._toggle_binding_residue("A:PHE:2A")
    composer = app._create_selected_variants()
    assert composer.setup == setup
    composer.add_button.invoke()
    assert affinity.variants.get("1.0", "end-1c").splitlines() == [
        "A_F2AA\tA:F2AA", "A_F2AA_2\tA:F2AA",
    ]
    assert not errors


@pytest.mark.parametrize("provider", ["affinity", "hotspot"])
def test_partner_two_surface_keeps_original_groups_when_appending_to_a_draft(desktop, tmp_path, provider):
    app, errors = desktop
    source, _multichain = _synthetic_sources(tmp_path)
    setup = {"structure": str(source), "chains_a": "A", "chains_b": "B"}
    affinity = app.open_binding_task("mutation", setup)
    affinity.window.withdraw()
    affinity.variants.insert("1.0", "existing\tA:F2AA")
    if provider == "affinity":
        analyze(str(source), "A", "B", output=str(tmp_path / "mutation"), task="mutation",
                mutations=["A:F2AA"], device="cpu", log=lambda _: None)
        affinity.load(tmp_path / "mutation" / "analysis.json")
        app.affinity_partner.current(1)
        app.affinity_partner.event_generate("<<ComboboxSelected>>")
        app.root.update()
    else:
        ranking = run_ranking(RankingSettings(str(source), "A", "B", "B", "", str(tmp_path / "ranking")))
        hotspot = app.open_binding_task("hotspot", setup)
        hotspot.window.withdraw()
        hotspot.display(ranking)
        hotspot.table.selection_set("B:PHE:2")
        assert hotspot.show_surface()
    assert app._successful_single_run["params"]["chain_a"] == "B"
    app._toggle_binding_residue("B:PHE:2")
    composer = app._create_selected_variants()
    assert composer.setup == setup
    composer.add_button.invoke()
    assert affinity.variants.get("1.0", "end-1c").splitlines() == [
        "existing\tA:F2AA", "B_F2A\tB:F2A",
    ]
    assert not errors


@pytest.mark.parametrize("side", [0, 1])
@pytest.mark.parametrize("multichain", [False, True])
def test_plugin_generated_atlas_restores_structure_for_followup_analysis(
    desktop, tmp_path, monkeypatch, side, multichain,
):
    app, errors = desktop
    if multichain:
        _single, source = _synthetic_sources(tmp_path)
        groups = [["A", "B"], ["C"]]
        selected_key = "A:PHE:2" if side == 0 else "C:PHE:2"
        mutation = "A:F2A" if side == 0 else "C:F2A"
    else:
        source = tmp_path / "structure.pdb"
        original = (Path(__file__).parent / "fixtures" / "tiny_complex.pdb").read_text()
        source.write_text("".join(
            line[:26] + "A" + line[27:]
            if line.startswith("ATOM  ") and line[21] == "A" and int(line[22:26]) == 1 else line
            for line in original.splitlines(keepends=True)
        ))
        groups = [["A"], ["B"]]
        selected_key = "A:GLY:1A" if side == 0 else "B:GLY:1"
        mutation = "A:G1AA" if side == 0 else "B:G1A"
    interface = topoppi.get_interface(topoppi.load_structure(source), *groups)
    coords_a, atoms_a, coords_b, atoms_b = interface.chain_data(side)
    patches, original_visualizer = make_atlas()
    visualizer = InterfaceVisualizer(
        atoms_a, coords_a, coords_b, atoms_b, chain_a_id=",".join(groups[side]),
        chain_b_id=",".join(groups[1 - side]), config=original_visualizer.config,
    )
    atlas = tmp_path / "partner.atlas.npz"
    save_atlas(atlas, patches, visualizer, run_metadata={
        "analysis_id": "synthetic-analysis", "side": side,
        "interface": {"source": {"path": str(source)}, "groups": groups},
    })
    monkeypatch.setattr("tkinter.filedialog.askopenfilename", lambda **_: str(atlas))
    app.open_atlas()
    assert app._successful_single_run["params"]["chain_a"] == ",".join(groups[side])
    assert app.binding_setup(displayed=True) == {
        "structure": str(source), "chains_a": ",".join(groups[0]), "chains_b": ",".join(groups[1]),
    }
    app._toggle_binding_residue(selected_key)
    composer = app._create_selected_variants()
    assert str(composer.add_button.cget("state")) == "normal"
    assert composer.preview.item(composer.preview.get_children()[0], "values")[1] == mutation
    composer.window.destroy()
    assert not errors


@pytest.mark.parametrize("provider", ["affinity", "hotspot"])
def test_saved_plugin_snapshot_recovers_original_groups_for_an_existing_draft(
    desktop, tmp_path, monkeypatch, provider,
):
    app, errors = desktop
    source, _multichain = _synthetic_sources(tmp_path)
    setup = {"structure": str(source), "chains_a": "A", "chains_b": "B"}
    affinity = app.open_binding_task("mutation", setup)
    affinity.window.withdraw()
    affinity.variants.insert("1.0", "existing\tA:F2AA")

    interface = topoppi.get_interface(topoppi.load_structure(source), "A", "B")
    coords_a, atoms_a, coords_b, atoms_b = interface.chain_data(1)
    patches, original_visualizer = make_atlas()
    visualizer = InterfaceVisualizer(
        atoms_a, coords_a, coords_b, atoms_b, chain_a_id="B", chain_b_id="A",
        config=original_visualizer.config,
    )
    # Earlier Save Atlas documents stored source groups in the plugin snapshot.
    # The numerical values here are arbitrary synthetic display data.
    if provider == "affinity":
        snapshot = {
            "analysis_json": str(tmp_path / "analysis" / "analysis.json"),
            "analysis_id": "synthetic-analysis", "input": {"path": str(source)},
            "model": {"target": "ddg"}, "groups": [["A"], ["B"]],
            "predictions": [{"variant_id": "existing", "mutations": "A:F2AA", "predicted": .25, "status": "ok"}],
            "variant_id": "existing", "side": 1,
        }
    else:
        snapshot = {
            "directory": str(tmp_path / "ranking"),
            "metadata": {"structure": str(source), "partner_groups": [["A"], ["B"]], "selected_partner": ["B"]},
            "rows": [{"residue_key": "B:PHE:2", "chain": "B", "residue_number": 2,
                      "insertion_code": "", "wt": "F", "rank_score": .25, "rank_min": 1, "rank_max": 1}],
            "selected": "B:PHE:2", "selected_residues": ["B:PHE:2"],
        }
    atlas = tmp_path / "saved.atlas.npz"
    save_atlas(atlas, patches, visualizer, run_metadata={
        "params": {"path": str(source), "chain_a": "B", "chain_b": "A", "min_points": 1},
        "manifest": {"config": {}, provider: snapshot},
    })
    monkeypatch.setattr("tkinter.filedialog.askopenfilename", lambda **_: str(atlas))
    app.open_atlas()
    assert app._successful_single_run["params"]["chain_a"] == "B"
    app._toggle_binding_residue("B:PHE:2")
    composer = app._create_selected_variants()
    assert composer.setup == setup
    composer.add_button.invoke()
    assert affinity.variants.get("1.0", "end-1c").splitlines() == [
        "existing\tA:F2AA", "B_F2A\tB:F2A",
    ]
    assert not errors
