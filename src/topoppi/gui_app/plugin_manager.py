"""Discover installed analysis tools and install optional plugin packages."""
from __future__ import annotations

import importlib
import json
import queue
import subprocess
import sys
import threading
import tkinter as tk
import webbrowser
from importlib.resources import files
from pathlib import Path
from tkinter import filedialog, ttk

from topoppi import __version__, plugins


def plugin_catalog():
    """Read display/install suggestions; runtime discovery remains entry-point based."""
    return {item["name"]: item for item in json.loads(files("topoppi").joinpath("assets/plugins.json").read_text())}


def installation_command(requirement):
    """Use the running application environment and retain its Core version."""
    return [sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
            "--progress-bar", "off", requirement, f"topoppi=={__version__}"]


class PluginManager:
    def __init__(self, host):
        self.host = host
        self.window = tk.Toplevel(host.root)
        self.window.title("TopoPPI Plugins")
        self.window.geometry("720x510")
        self.window.minsize(590, 430)
        self.window.transient(host.root)
        self.window.protocol("WM_DELETE_WINDOW", self.close)
        self.events = queue.Queue()
        self.process = None
        self.working = False
        self.cancelled = False
        self.catalog = plugin_catalog()
        self.rows = {}
        body = ttk.Frame(self.window, padding=18)
        body.pack(fill=tk.BOTH, expand=True)
        ttk.Label(body, text="Tools for your interface", font=("TkDefaultFont", 16, "bold")).pack(anchor=tk.W)
        ttk.Label(body, text="Choose an analysis, open its workspace, or add a plugin.",
                  wraplength=650).pack(anchor=tk.W, pady=(5, 14))
        self.table = ttk.Treeview(body, columns=("name", "state", "version"), show="headings", height=5)
        for key, title, width in (("name", "Plugin", 260), ("state", "Status", 150), ("version", "Version", 120)):
            self.table.heading(key, text=title)
            self.table.column(key, width=width, minwidth=80)
        self.table.pack(fill=tk.BOTH, expand=True)
        self.table.bind("<<TreeviewSelect>>", self._select)
        self.table.bind("<Double-1>", lambda _event: self.open_selected())
        self.description = tk.StringVar()
        self.detail_label = ttk.Label(body, textvariable=self.description, wraplength=650, justify=tk.LEFT)
        self.detail_label.pack(fill=tk.X, pady=(12, 8))
        actions = ttk.Frame(body)
        actions.pack(fill=tk.X)
        self.open_button = ttk.Button(actions, text="Open", command=self.open_selected)
        self.open_button.pack(side=tk.LEFT)
        self.install_button = ttk.Button(actions, text="Install", command=self.install_selected)
        self.install_button.pack(side=tk.LEFT, padx=6)
        self.guide_button = ttk.Button(actions, text="User guide", command=self.open_guide)
        self.guide_button.pack(side=tk.LEFT)
        self.refresh_button = ttk.Button(actions, text="Refresh", command=self.refresh)
        self.refresh_button.pack(side=tk.RIGHT)
        bottom = ttk.Frame(body)
        bottom.pack(fill=tk.X, pady=(10, 0))
        self.local_button = ttk.Button(bottom, text="Install from file...", command=self.install_file)
        self.local_button.pack(side=tk.LEFT)
        self.cancel_button = ttk.Button(bottom, text="Cancel installation", command=self.cancel, state=tk.DISABLED)
        self.status = tk.StringVar(value="Chemical interactions are included. Add prediction tools when you need them.")
        self.status_label = ttk.Label(body, textvariable=self.status, wraplength=650, justify=tk.LEFT)
        self.status_label.pack(fill=tk.X, pady=(12, 3))
        self.progress = ttk.Progressbar(body, mode="indeterminate")
        self.show_log = tk.BooleanVar(value=False)
        self.log_toggle = ttk.Checkbutton(body, text="Installation details", variable=self.show_log,
                                         command=self._toggle_log)
        self.log = tk.Text(body, height=6, wrap=tk.WORD, state=tk.DISABLED)
        body.bind("<Configure>", lambda event: [label.configure(wraplength=max(200, event.width))
                                                for label in (self.detail_label, self.status_label)])
        self.refresh(rediscover=False)
        self.poll_id = self.window.after(100, self._poll)

    def _toggle_log(self):
        if self.show_log.get():
            self.log.pack(fill=tk.BOTH, expand=True, pady=(4, 0))
        else:
            self.log.pack_forget()

    def refresh(self, *, rediscover=True):
        if self.working:
            return
        if rediscover:
            importlib.invalidate_caches()
            plugins.refresh()
            self.host._init_extensions()
        selected = self.table.selection()
        self.rows = {key: dict(value, installed=False, version="") for key, value in self.catalog.items()}
        for info in plugins.list():
            row = self.rows.setdefault(info.name, {"name": info.name, "title": info.name.replace("_", " ").title()})
            try:
                info = plugins.info(info.name)
                row.update(installed=True, version=info.version)
                row.setdefault("description", info.description)
            except Exception as exc:
                row.update(installed=True, version=info.version, error=str(exc))
        self.table.delete(*self.table.get_children())
        errors = getattr(self.host, "_extension_errors", {})
        for name, row in self.rows.items():
            if name in errors:
                row["error"] = errors[name]
            state = ("Setup needed" if row.get("error") else
                     "Included" if row.get("installed") and row.get("included") else
                     "Installed" if row.get("installed") else "Available to install")
            self.table.insert("", tk.END, iid=name, values=(row["title"], state, row["version"]))
        if self.rows:
            self.table.selection_set(selected[0] if selected and selected[0] in self.rows else next(iter(self.rows)))
        self._select()

    def _selected(self):
        selected = self.table.selection()
        return self.rows.get(selected[0], {}) if selected else {}

    def _select(self, _event=None):
        row = self._selected()
        detail = row.get("description", "")
        if row.get("error"):
            detail += "\nSetup details: " + row["error"]
        available = row.get("name") in self.host._gui_extensions
        if row.get("installed") and not available and not row.get("error"):
            detail += "\nUse this plugin through its command-line or Python interface."
        self.description.set(detail)
        self.open_button.configure(state=tk.NORMAL if available and not self.working and not self.host._busy
                                   else tk.DISABLED)
        self.install_button.configure(state=tk.NORMAL if row.get("install") and not row.get("installed")
                                      and not self.working else tk.DISABLED)
        self.guide_button.configure(state=tk.NORMAL if row.get("guide") else tk.DISABLED)

    def open_selected(self):
        row = self._selected()
        if self.working or self.host._busy or row.get("name") not in self.host._gui_extensions:
            return
        try:
            self.host.open_extension(row["name"])
            self.window.withdraw()
        except Exception as exc:
            self.status.set(f"Open {row['title']}: {exc}")

    def open_guide(self):
        if url := self._selected().get("guide"):
            webbrowser.open(url)

    def install_selected(self):
        row = self._selected()
        if row.get("install") and not row.get("installed"):
            self.install(row["install"], row["title"])

    def install_file(self):
        path = filedialog.askopenfilename(parent=self.window, title="Choose a plugin package",
                                          filetypes=[("Python wheel", "*.whl"), ("Source archive", "*.tar.gz")])
        if path:
            self.install(path, Path(path).name)

    def install(self, requirement, title):
        if self.working or self.host._busy:
            self.status.set("Finish the current task, then install the plugin.")
            return
        self.working = True
        self.cancelled = False
        self.status.set(f"Installing {title} in this TopoPPI environment...")
        self.progress.pack(fill=tk.X, after=self.status_label)
        self.progress.start(15)
        self.local_button.configure(state=tk.DISABLED)
        self.refresh_button.configure(state=tk.DISABLED)
        self.cancel_button.configure(state=tk.NORMAL)
        self.cancel_button.pack(side=tk.RIGHT)
        self.log_toggle.pack(anchor=tk.W, pady=(6, 0))
        self._select()
        threading.Thread(target=self._install_worker, args=(requirement,), daemon=True).start()

    def _install_worker(self, requirement):
        try:
            options = {"creationflags": subprocess.CREATE_NO_WINDOW} if sys.platform == "win32" else {}
            self.process = subprocess.Popen(installation_command(requirement), stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                                            **options)
            if self.cancelled:
                self.process.terminate()
            for line in self.process.stdout:
                self.events.put(("line", line))
            self.events.put(("done", self.process.wait()))
        except OSError as exc:
            self.events.put(("line", str(exc)))
            self.events.put(("done", 1))

    def cancel(self):
        self.cancelled = True
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
        self.status.set("Cancelling installation...")

    def _poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "line":
                    self.log.configure(state=tk.NORMAL)
                    self.log.insert(tk.END, value)
                    self.log.see(tk.END)
                    self.log.configure(state=tk.DISABLED)
                else:
                    self.working = False
                    self.progress.stop()
                    self.progress.pack_forget()
                    self.cancel_button.configure(state=tk.DISABLED)
                    self.cancel_button.pack_forget()
                    self.local_button.configure(state=tk.NORMAL)
                    self.refresh_button.configure(state=tk.NORMAL)
                    self.refresh()
                    self.status.set("Installation cancelled. Refresh to check the available tools." if self.cancelled else
                                    "Plugin installed. Select it and choose Open." if value == 0 else
                                    "Installation needs attention. Open Installation details or choose a local package.")
        except queue.Empty:
            pass
        self.poll_id = self.window.after(100, self._poll)

    def close(self):
        if self.working:
            self.status.set("Installation is running. Choose Cancel installation to stop it.")
            return False
        self.window.after_cancel(self.poll_id)
        self.window.destroy()
        return True
