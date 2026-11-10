"""Small optional GUI extension lifecycle and task-neutral saved metadata."""
from __future__ import annotations

import json
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk

from topoppi import __version__

from .binding_analysis import BindingAnalysisMixin


class ExtensionsMixin(BindingAnalysisMixin):
    def _init_extensions(self):
        from topoppi import plugins

        if not hasattr(self, "_extension_errors"):
            self._extension_errors = {}
        if hasattr(self, "plugin_actions_frame"):
            self._init_binding_analysis()
        if hasattr(self, "plugin_actions_frame") and not hasattr(self, "btn_plugins"):
            self.btn_plugins = ttk.Button(self.plugin_actions_frame, text="Plugins...", command=self.show_plugins)
            self.btn_plugins.pack(side=tk.LEFT, padx=(8, 0))
        for action in plugins.gui_actions():
            if action.name in self._gui_extensions:
                continue
            try:
                controller = plugins.open_gui(action.name, self)
                initialize = getattr(controller, "initialize", None)
                if initialize:
                    initialize()
                self._gui_extensions[action.name] = controller
                self._extension_errors.pop(action.name, None)
            except Exception as exc:
                self._extension_errors[action.name] = str(exc)
                self.log_history.append(f"Plugin interface {action.name} is unavailable: {exc}")

    def get_current_geometry(self):
        """Return native geometry associated with the displayed result, if retained."""
        run = self._successful_single_run
        return run.get("geometry") if run else None

    def show_plugins(self):
        from topoppi.gui_app.plugin_manager import PluginManager

        manager = getattr(self, "_plugin_manager", None)
        if manager is not None and manager.window.winfo_exists():
            manager.window.deiconify()
            manager.window.lift()
            manager.refresh(rediscover=False)
        else:
            self._plugin_manager = PluginManager(self)

    def open_extension(self, name, setup=None):
        controller = self._gui_extensions.get(name)
        if controller is None:
            raise RuntimeError(f"GUI extension '{name}' is unavailable. Install its corresponding TopoPPI plugin.")
        return controller.open(setup=setup)

    def _sync_extensions(self):
        binding_window = getattr(self, "_binding_window", None)
        if binding_window is not None and binding_window.winfo_exists():
            self._binding_launcher_refresh()
        manager = getattr(self, "_plugin_manager", None)
        if manager is not None and manager.window.winfo_exists():
            manager._select()
        for controller in self._gui_extensions.values():
            callback = getattr(controller, "sync_controls", None)
            if callback:
                callback()

    def _extensions_need_restart(self):
        return any(getattr(controller, "restart_required", False) for controller in self._gui_extensions.values())

    def _extensions_redraw(self, style):
        for controller in self._gui_extensions.values():
            callback = getattr(controller, "redraw", None)
            if callback and callback(style):
                return True
        return False

    def _extensions_document_loaded(self, manifest):
        # Preserve unknown metadata as plain data: viewing it imports no plugin.
        self._saved_extension_metadata = manifest
        for controller in self._gui_extensions.values():
            callback = getattr(controller, "document_loaded", None)
            if callback:
                callback(manifest)

    def _extension_filename(self, run):
        for controller in self._gui_extensions.values():
            callback = getattr(controller, "filename", None)
            if callback and (name := callback(run)):
                return name
        return None

    def _close_extensions(self):
        manager = getattr(self, "_plugin_manager", None)
        if manager is not None and manager.window.winfo_exists() and not manager.close():
            return False
        for controller in self._gui_extensions.values():
            callback = getattr(controller, "close", None)
            if callback and not callback():
                return False
        return True

    def saved_metadata_text(self):
        run = self._successful_single_run or self._pending_single_run
        metadata = run["manifest"] if run else self._saved_extension_metadata
        return json.dumps(metadata, indent=2, default=str)

    def show_saved_metadata(self):
        window = tk.Toplevel(self.root)
        window.title("Saved analysis metadata")
        text = scrolledtext.ScrolledText(window, width=88, height=28, wrap=tk.WORD)
        text.pack(fill=tk.BOTH, expand=True)
        text.insert("1.0", self.saved_metadata_text())
        text.configure(state=tk.DISABLED)

    def diagnostics_text(self):
        run = self._successful_single_run or self._pending_single_run
        extensions = {}
        for name, controller in self._gui_extensions.items():
            callback = getattr(controller, "diagnostics", None)
            if callback:
                extensions[name] = callback()
        details = {"topoppi": __version__, "python": sys.version, "executable": sys.executable,
                   "input": {"path": self.var_input_path.get(), "chain_a": self.var_chain_a.get(),
                             "chain_b": self.var_chain_b.get()},
                   "displayed_run": run["manifest"] if run else None, "extensions": extensions}
        return json.dumps(details, indent=2, default=str) + "\n\nActivity\n" + "\n".join(self.log_history)

    def copy_diagnostics(self):
        self.root.clipboard_clear()
        self.root.clipboard_append(self.diagnostics_text())
        self.status_var.set("Run settings and activity copied to the clipboard.")

    def save_diagnostics(self):
        path = filedialog.asksaveasfilename(parent=self.root, title="Save diagnostics", initialfile="topoppi-diagnostics.txt",
                                           defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if path:
            try:
                Path(path).write_text(self.diagnostics_text(), encoding="utf-8")
                self.status_var.set(f"Diagnostics saved to {path}")
            except OSError as exc:
                messagebox.showerror("Save diagnostics", str(exc), parent=self.root)
