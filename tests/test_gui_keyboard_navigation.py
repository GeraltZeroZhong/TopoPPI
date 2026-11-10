"""Keyboard access to controls beyond the visible main-window sidebar."""

import os
import sys
import tkinter as tk
from unittest.mock import patch

import pytest

from topoppi.gui_app.app import InterfaceApp

pytestmark = pytest.mark.skipif(
    sys.platform not in {"win32", "darwin"} and not os.environ.get("DISPLAY"),
    reason="A graphical display is required",
)


@pytest.fixture
def desktop(tmp_path):
    with patch.object(InterfaceApp, "_recent_store_path", return_value=tmp_path / "recent.json"):
        root = tk.Tk()
        app = InterfaceApp(root)
        root.geometry("960x600")
        root.update()
        errors = []
        root.report_callback_exception = lambda kind, value, traceback: errors.append(value)
        try:
            yield app
            assert not errors
        finally:
            app.close()


@pytest.mark.parametrize(
    ("page", "first_name", "target_name"),
    [
        ("basic", "entry_file", "entry_output_dir_basic"),
        ("advanced", "entry_file_advanced", "entry_cutoff"),
    ],
)
def test_tab_reveals_sidebar_control_and_preserves_input_focus(desktop, page, first_name, target_name):
    app = desktop
    app.var_settings_page.set(page)
    app._sync_settings_page()
    app.sidebar_canvas.yview_moveto(0)
    first, target = getattr(app, first_name), getattr(app, target_name)
    first.focus_force()
    app.root.update()

    for _ in range(40):
        focused = app.root.focus_get()
        if focused is target:
            break
        focused.event_generate("<Tab>")
        app.root.update()

    assert app.root.focus_get() is target
    viewport_top = app.sidebar_canvas.winfo_rooty()
    viewport_bottom = viewport_top + app.sidebar_canvas.winfo_height()
    assert target.winfo_rooty() >= viewport_top
    assert target.winfo_rooty() + target.winfo_height() <= viewport_bottom


def test_page_keys_scroll_sidebar_while_entry_keeps_focus(desktop):
    app = desktop
    app.var_settings_page.set("advanced")
    app._sync_settings_page()
    app.sidebar_canvas.yview_moveto(0)
    app.entry_file_advanced.focus_force()
    app.root.update()
    before = app.sidebar_canvas.yview()[0]

    app.entry_file_advanced.event_generate("<Next>")
    app.root.update()
    after = app.sidebar_canvas.yview()[0]
    assert after > before
    assert app.root.focus_get() is app.entry_file_advanced

    app.entry_file_advanced.event_generate("<Prior>")
    app.root.update()
    assert app.sidebar_canvas.yview()[0] < after
    assert app.root.focus_get() is app.entry_file_advanced
