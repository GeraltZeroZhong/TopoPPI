"""Scrolling settings must not silently change the next analysis or its style."""

import os
import sys
import tkinter as tk
from pathlib import Path
from unittest.mock import patch

import pytest

from topoppi import plugins
from topoppi.gui_app.app import InterfaceApp

pytestmark = pytest.mark.skipif(
    sys.platform not in {"win32", "darwin"} and not os.environ.get("DISPLAY"),
    reason="A graphical display is required",
)


@pytest.fixture
def desktop(tmp_path, monkeypatch):
    monkeypatch.setattr(plugins, "_default_registry", plugins.PluginRegistry(entry_points=()))
    with patch.object(InterfaceApp, "_recent_store_path", return_value=tmp_path / "recent.json"):
        root = tk.Tk()
        app = InterfaceApp(root)
        root.geometry("960x600")
        root.update()
        errors = []
        root.report_callback_exception = lambda _kind, value, _traceback: errors.append(value)
        try:
            yield app
            assert not errors
        finally:
            app.close()


@pytest.mark.parametrize("sequence,options", [
    ("<MouseWheel>", {"delta": -120}),
    ("<Button-5>", {}),
])
def test_scrolling_over_chain_field_keeps_selected_surface(desktop, sequence, options):
    app = desktop
    app.var_input_path.set(str(Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"))
    app._commit_structure_path()
    app.root.update()
    app.sidebar_canvas.yview_moveto(0)
    app.sidebar_canvas.event_generate("<Enter>")
    before = app.sidebar_canvas.yview()[0]

    app.combo_basic_chain_a.event_generate(sequence, **options)
    app.root.update()

    assert app.sidebar_canvas.yview()[0] > before
    assert app.var_chain_a.get() == "A"
    assert app.var_chain_b.get() == "B"
    assert "disabled" not in app.btn_run.state()


def test_scrolling_over_font_size_keeps_numeric_setting(desktop):
    app = desktop
    app.var_settings_page.set("appearance")
    app._sync_settings_page()
    app.spin_size.focus_force()
    app.root.update()
    before_scroll = app.sidebar_canvas.yview()[0]
    before_size = app.spin_size.get()

    app.spin_size.event_generate("<MouseWheel>", delta=120)
    app.root.update()

    assert app.sidebar_canvas.yview()[0] < before_scroll
    assert app.spin_size.get() == before_size


@pytest.mark.skipif(sys.platform == "darwin", reason="macOS uses a native combobox menu")
def test_keyboard_can_still_choose_a_chain(desktop):
    app = desktop
    combo = app.combo_basic_chain_a
    combo.configure(values=("A", "B", "C"))
    combo.set("A")
    combo.focus_force()
    app.root.update()
    combo.event_generate("<Down>")
    app.root.update()
    popup = app.root.tk.call("ttk::combobox::PopdownWindow", str(combo))
    listbox = popup + ".f.l"

    app.root.tk.call("event", "generate", listbox, "<Down>")
    app.root.tk.call("event", "generate", listbox, "<Return>")
    app.root.update()

    assert combo.get() == "B"
    assert not app.root.tk.call("winfo", "ismapped", popup)


@pytest.mark.skipif(sys.platform == "darwin", reason="macOS uses a native combobox menu")
def test_open_chain_list_retains_its_own_scroll(desktop):
    app = desktop
    combo = app.combo_basic_chain_a
    combo.configure(values=tuple(f"C{index}" for index in range(40)))
    combo.set("C0")
    combo.focus_force()
    app.root.update()
    combo.event_generate("<Down>")
    app.root.update()
    popup = app.root.tk.call("ttk::combobox::PopdownWindow", str(combo))
    listbox = popup + ".f.l"
    before_list = app.root.tk.call(listbox, "yview")[0]
    before_sidebar = app.sidebar_canvas.yview()

    try:
        app.root.tk.call("event", "generate", listbox, "<MouseWheel>", "-delta", -120)
        app.root.update()
        assert app.root.tk.call(listbox, "yview")[0] > before_list
        assert app.sidebar_canvas.yview() == before_sidebar
        assert combo.get() == "C0"
    finally:
        app.root.tk.call("event", "generate", listbox, "<Escape>")
        app.root.update()
