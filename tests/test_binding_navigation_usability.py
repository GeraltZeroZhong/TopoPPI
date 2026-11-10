"""Binding task navigation follows installed tools and preserves actual drafts."""
import tkinter as tk
from pathlib import Path

import pytest

from topoppi.gui_app.binding_analysis import BindingAnalysisMixin


@pytest.fixture
def root():
    try:
        window = tk.Tk()
    except tk.TclError:
        pytest.skip("A graphical display is required")
    window.withdraw()
    yield window
    window.destroy()


def _labels(parent):
    result = []
    for child in parent.winfo_children():
        if child.winfo_class() == "TLabel":
            result.append(str(child.cget("text")))
        result.extend(_labels(child))
    return result


def test_reopening_launcher_refreshes_plugin_install_guidance(root):
    host = BindingAnalysisMixin()
    host.root = root
    host._binding_tasks = {}
    host._successful_single_run = None
    host.var_input_path = tk.StringVar(root, "")
    host.var_chain_a = tk.StringVar(root, "A")
    host.var_chain_b = tk.StringVar(root, "B")
    host.show_plugins = lambda: None
    host.show_binding_analysis()
    root.update()
    assert any("Install the Hotspot plugin" in text for text in _labels(host._binding_window))

    # PluginManager.refresh registers newly available tasks on the same host.
    host.register_binding_task("hotspot", title="Rank hotspots", description="", plugin="hotspot", open_task=lambda **_: None)
    host.show_binding_analysis()
    root.update()
    assert not any("Install the Hotspot plugin" in text for text in _labels(host._binding_window))
    assert any("Install the Affinity plugin" in text for text in _labels(host._binding_window))


def test_switching_structure_without_variants_does_not_report_a_nonexistent_draft(root, tmp_path):
    gui = pytest.importorskip("topoppi_affinity.gui")
    panel = gui.AffinityPanel(root)
    source = Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"
    other = tmp_path / "next-complex.pdb"
    other.write_text(source.read_text())
    try:
        panel.vars["structure"].set(str(source))
        panel.vars["chains_a"].set("A")
        panel.vars["chains_b"].set("B")
        panel.restore_setup({"prediction_task": "absolute", "structure": str(other), "chains_a": "A", "chains_b": "B"})
        root.update()
        assert panel.vars["structure"].get() == str(other)
        assert panel.task.get() == "Binding affinity (pKd)"
        assert panel.worker is None

        # Actual mutation drafts must still prevent silently retargeting residues.
        panel.variants.insert("1.0", "candidate\tA:A1V")
        with pytest.raises(ValueError, match="different structure"):
            panel.restore_setup({"prediction_task": "mutation", "structure": str(source), "chains_a": "A", "chains_b": "B"})
        assert panel.vars["structure"].get() == str(other)
        assert panel.variants.get("1.0", "end-1c") == "candidate\tA:A1V"
    finally:
        panel.close()
