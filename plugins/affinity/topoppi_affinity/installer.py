"""Install the optional affinity plugin into the running TopoPPI environment.

This module uses only the standard library. Importing or opening it never starts
pip or downloads models; the CLI invocation or the Install button does that.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import importlib.util
import json
import os
import queue
import re
import signal
import subprocess
import sys
import tempfile
import threading
from collections import deque
from pathlib import Path

from topoppi._version import __version__

TASKS = ("absolute", "mutation", "known", "unseen")
TASK_LABELS = {
    "absolute": "Binding affinity (pKd)",
    "mutation": "Mutation effects (bundled CPU model)",
    "known": "Mutation effects: known system",
    "unseen": "Mutation effects: unseen system",
}


class InstallationCancelled(RuntimeError):
    pass


def installation_status() -> dict[str, str | None]:
    try:
        plugin = importlib.metadata.version("topoppi-affinity")
    except importlib.metadata.PackageNotFoundError:
        plugin = None
    return {"python": sys.executable, "core": __version__, "plugin": plugin}


def _plugin_directory(source: str | Path) -> Path:
    folder = Path(source).expanduser().resolve()
    if (folder / "plugins" / "affinity" / "pyproject.toml").is_file():
        folder = folder / "plugins" / "affinity"
    if not (folder / "pyproject.toml").is_file():
        raise ValueError(f"No plugin pyproject.toml in {folder}; choose plugins/affinity or its TopoPPI checkout.")
    metadata = (folder / "pyproject.toml").read_text(encoding="utf-8")
    if not re.search(r"^name\s*=\s*[\"']topoppi-affinity[\"']\s*$", metadata, re.MULTILINE):
        raise ValueError(f"{folder} is not the topoppi-affinity project; choose its plugins/affinity directory.")
    return folder


def discover_plugin_source() -> Path | None:
    """Find a plugin in this source checkout or the current working directory."""
    for start in (Path(__file__).resolve().parent, Path.cwd()):
        for candidate in (start, *start.parents):
            if not ((candidate / "plugins" / "affinity" / "pyproject.toml").is_file()
                    or (candidate / "pyproject.toml").is_file()):
                continue
            try:
                return _plugin_directory(candidate)
            except ValueError:
                continue
    return None


def plugin_requirement(source: str | Path | None = None, *, models=False) -> str:
    extra = "[models]" if models else ""
    folder = _plugin_directory(source) if source is not None else discover_plugin_source()
    if folder is not None:
        return f"{folder}{extra}"
    if ".dev" in __version__:
        raise ValueError("This development version needs the matching local plugin source. "
                         "Choose your TopoPPI checkout in Locations, or pass --source /path/to/TopoPPI. "
                         "A published package for this development version is not assumed to exist.")
    return f"topoppi-affinity{extra}=={__version__}"


def _selected_tasks(tasks=None) -> tuple[str, ...]:
    selected = tuple(dict.fromkeys(TASKS[:1] if tasks is None else tasks))
    if not selected or any(task not in TASKS for task in selected):
        raise ValueError("Choose a prediction task: absolute, mutation, known or unseen.")
    return selected


def _stop_process(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    if os.name == "nt":
        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
                       creationflags=subprocess.CREATE_NO_WINDOW)
    else:
        os.killpg(process.pid, signal.SIGTERM)
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            process.kill()
        else:
            os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def _run(command: list[str], progress, cancel_event: threading.Event) -> None:
    if cancel_event.is_set():
        raise InstallationCancelled("Preparation cancelled. You can retry; completed downloads remain cached.")
    options = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {"start_new_session": True}
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, encoding="utf-8", errors="replace", bufsize=1,
                               env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"}, **options)
    output = queue.Queue()
    tail = deque(maxlen=14)

    def read_output():
        try:
            for line in process.stdout:
                output.put(line.rstrip())
        finally:
            output.put(None)

    reader = threading.Thread(target=read_output, daemon=True)
    reader.start()
    try:
        while True:
            if cancel_event.is_set():
                _stop_process(process)
                raise InstallationCancelled("Preparation cancelled. Retry to finish installation or downloads.")
            try:
                line = output.get(timeout=0.1)
            except queue.Empty:
                continue
            if line is None:
                break
            tail.append(line)
            progress(line)
        while process.poll() is None:
            if cancel_event.wait(0.1):
                _stop_process(process)
                raise InstallationCancelled("Preparation cancelled. Retry to finish installation or downloads.")
        returncode = process.returncode
        if returncode:
            raise RuntimeError(f"Preparation command failed (exit {returncode}):\n" + "\n".join(tail))
    except BaseException:
        _stop_process(process)
        raise
    finally:
        reader.join(timeout=2)
        process.stdout.close()


def install_affinity(*, source=None, download=True, tasks=None, cache_dir=None, progress=None,
                     cancel_event=None, install_packages=True, on_install_started=None, on_installed=None) -> None:
    """Install compatible dependencies and prepare resources for selected tasks.

    The default task is absolute affinity. The running core version is retained;
    resource preparation can be retried independently with install_packages=False.
    """
    progress = progress or (lambda text: print(text, file=sys.stderr, flush=True))
    cancel_event = cancel_event or threading.Event()
    selected = _selected_tasks(tasks) if download or tasks else ()
    if not install_packages and not download:
        raise ValueError("Choose package installation or model preparation.")
    requirement = (plugin_requirement(source, models=True) if any(t in {"known", "unseen"} for t in selected)
                   else plugin_requirement(source)) if install_packages else None
    packages_ready = not install_packages
    progress(f"TopoPPI {__version__}; Python: {sys.executable}")
    try:
        if install_packages:
            if cancel_event.is_set():
                raise InstallationCancelled("Preparation cancelled. Retry when ready.")
            progress(f"Installing {requirement}; the existing TopoPPI version will be retained.")
            if on_install_started:
                on_install_started()
            if importlib.util.find_spec("pip") is None:
                _run([sys.executable, "-m", "ensurepip"], progress, cancel_event)
            with tempfile.TemporaryDirectory(prefix="topoppi-affinity-install-") as directory:
                constraint = Path(directory) / "core-version.txt"
                constraint.write_text(f"topoppi=={__version__}\n", encoding="utf-8")
                _run([sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
                      "--constraint", str(constraint), requirement], progress, cancel_event)
            packages_ready = True
            if on_installed:
                on_installed()
            importlib.invalidate_caches()
        if download:
            progress("Preparing selected tasks: " + ", ".join(TASK_LABELS[task] for task in selected) + ".")
            if any(task in {"known", "unseen"} for task in selected):
                progress("Checking advanced model resources and reusing cached files. "
                         "ESMC sequence weights and neural dependencies can require several GB.")
            else:
                progress("Checking bundled CPU models; their weights are included with the plugin.")
            command = [sys.executable, "-m", "topoppi_affinity", "models", "download"]
            for task in selected:
                command.extend(["--task", task])
            if cache_dir:
                command.extend(["--cache-dir", str(Path(cache_dir).expanduser())])
            _run(command, progress, cancel_event)
    except InstallationCancelled:
        raise
    except (OSError, RuntimeError) as exc:
        retry = ("Retry model preparation with --models-only; completed downloads remain cached."
                 if packages_ready else "Installed packages and completed downloads are retained for retry.")
        raise RuntimeError(f"{exc}\nReview the reported cause before retrying. {retry}") from exc
    progress("Selected prediction resources are ready." if download else
             "Affinity plugin installed. Prepare prediction resources when needed.")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="topoppi-install-affinity",
                                     description="Install optional affinity prediction; prepare only the tasks you need.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--source", type=Path, help="Plugin directory or TopoPPI checkout; a local checkout is detected automatically.")
    tasks = parser.add_mutually_exclusive_group()
    tasks.add_argument("--task", action="append", choices=TASKS, help="Task to prepare; repeat for several tasks (default: absolute).")
    tasks.add_argument("--all", action="store_true", help="Prepare both bundled tasks: affinity and mutation effects.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--skip-download", action="store_true", help="Install packages now; prepare models later.")
    mode.add_argument("--models-only", action="store_true", help="Prepare resources using the installed plugin, without running pip.")
    parser.add_argument("--cache-dir", type=Path, help="Shared affinity model cache directory.")
    parser.add_argument("--status", action="store_true", help="Show versions and the target interpreter without changing anything.")
    args = parser.parse_args(argv)
    if args.status:
        for key, value in installation_status().items():
            print(f"{key}: {value or 'not installed'}")
        source = discover_plugin_source()
        print(f"local source: {source or 'not detected'}")
        return 0
    try:
        install_affinity(source=args.source, download=not args.skip_download, tasks=TASKS[:2] if args.all else args.task,
                         install_packages=not args.models_only, cache_dir=args.cache_dir)
    except KeyboardInterrupt:
        print("Installation interrupted. Retry to complete preparation.", file=sys.stderr)
        return 130
    except (ValueError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    return 0


def restart_command(host, *, task=None) -> list[str]:
    """Preserve the selected analysis input across a fresh-interpreter restart."""
    command = [sys.executable, "-m", "topoppi.gui", "--plugin", "affinity"]
    for field, plugin_field, flag in (("var_input_path", "structure", "--structure"),
                        ("var_chain_a", "chains_a", "-A"), ("var_chain_b", "chains_b", "-B"),
                        ("var_output_dir", "output", "--output-dir"), ("entry_optcuts_bin", "optcuts_bin", "--optcuts-bin")):
        widget = getattr(host, field, None) or getattr(host, "vars", {}).get(plugin_field)
        if widget is not None and widget.get().strip():
            command.extend([flag, widget.get().strip()])
    panel = host if hasattr(host, "export_setup") else getattr(host, "_affinity_panel", None)
    setup = panel.export_setup() if panel is not None and hasattr(panel, "export_setup") else None
    if setup is None and task is not None:
        setup = {"task": task}
    if setup is not None:
        with tempfile.NamedTemporaryFile(mode="w", prefix="topoppi-affinity-restart-", suffix=".json",
                                         encoding="utf-8", delete=False) as handle:
            json.dump(setup, handle)
            command.extend(["--plugin-setup", handle.name])
    return command


def open_installer(host, *, download=True, tasks=None):
    """Show task preparation without importing prediction dependencies or starting downloads."""
    import tkinter as tk
    from tkinter import filedialog, scrolledtext, ttk

    previous = getattr(host, "_affinity_installer", None)
    if previous is not None and previous.window.winfo_exists():
        previous.window.lift()
        return previous
    selected = _selected_tasks(tasks)

    class Installer:
        def __init__(self):
            self.window = tk.Toplevel(host.root)
            self.window.title("Prepare TopoPPI Affinity")
            self.window.geometry(f"760x{min(610, max(440, host.root.winfo_screenheight() - 80))}")
            self.window.minsize(500, 440)
            self.worker = None
            self.cancel = threading.Event()
            self.messages = queue.Queue()
            self.poll_id = None
            self.closing = False
            self.ready = bool(previous and previous.ready)
            owner = getattr(host, "host", None) or host
            self.restart_required = bool(getattr(owner, "_affinity_restart_required", False)
                                         or (previous and previous.restart_required))
            self.finished = True
            environment = installation_status()
            self.dependencies_installed = (previous.dependencies_installed if previous else
                                           bool(download and environment["plugin"] == environment["core"]
                                                and all(importlib.util.find_spec(module) is not None
                                                        for module in ("numpy", "pandas", "sklearn"))))
            detected = discover_plugin_source()
            self.source = tk.StringVar(value=previous.source.get() if previous else str(detected or ""))
            self.installed_source = previous.installed_source if previous else (self.source.get().strip() or None)
            self.cache = tk.StringVar(value=previous.cache.get() if previous else "")
            self.prepare_models = tk.BooleanVar(value=previous.prepare_models.get() if previous else download)
            self.task_vars = {task: tk.BooleanVar(value=(previous.task_vars[task].get() if previous and tasks is None
                                                         else task in selected)) for task in TASKS}
            self.inputs = []
            frame = ttk.Frame(self.window, padding=12)
            frame.pack(fill="both", expand=True)
            frame.columnconfigure(0, weight=1)
            frame.rowconfigure(1, weight=1)
            ttk.Label(frame, text="Choose the predictions you want to use", font=("TkDefaultFont", 12, "bold")).grid(
                row=0, column=0, sticky="w", pady=(0, 8))
            self.tabs = ttk.Notebook(frame)
            self.tabs.grid(row=1, column=0, sticky="nsew")
            choice = ttk.Frame(self.tabs, padding=10)
            locations = ttk.Frame(self.tabs, padding=10)
            output = ttk.Frame(self.tabs, padding=5)
            self.tabs.add(choice, text="Prediction tasks")
            self.tabs.add(locations, text="Locations")
            self.tabs.add(output, text="Progress")
            self.progress_tab = output
            models = ttk.Checkbutton(choice, text="Prepare prediction resources now", variable=self.prepare_models)
            models.pack(anchor="w", pady=(0, 5))
            self.inputs.append(models)
            for task in TASKS:
                widget = ttk.Checkbutton(choice, text=TASK_LABELS[task], variable=self.task_vars[task])
                widget.pack(anchor="w", padx=18, pady=2)
                self.inputs.append(widget)
            descriptions = [
                "Binding affinity and mutation effects include bundled CPU models, ready to use after installation.",
                "Advanced neural models use a separate model catalog. Known system: the interaction is represented in training. "
                "Unseen system: the interaction is outside the training systems.",
                "Advanced neural models can require several GB for ESMC sequence weights and prediction dependencies. "
                "Resource preparation reuses cached files.",
                "To install dependencies first, clear 'Prepare prediction resources now'. Save your figures before the required "
                "restart. Your structure, chains and prediction setup are retained.",
            ]
            for description in descriptions:
                ttk.Label(choice, text=description, wraplength=650).pack(fill="x", anchor="w", pady=(6, 0))
            ttk.Label(locations, text=f"TopoPPI {environment['core']} · Plugin: {environment['plugin'] or 'not installed'}").pack(anchor="w")
            ttk.Label(locations, text=f"Python environment: {environment['python']}", wraplength=650).pack(anchor="w", pady=(4, 8))
            for label, variable in (("Model storage folder (leave blank to use the shared cache)", self.cache),
                                    ("Plugin source (local checkout detected automatically)", self.source)):
                section = ttk.LabelFrame(locations, text=label, padding=6)
                section.pack(fill="x", pady=5)
                entry = ttk.Entry(section, textvariable=variable)
                entry.pack(side="left", fill="x", expand=True)
                browse = ttk.Button(section, text="Browse...", command=lambda var=variable: var.set(
                    filedialog.askdirectory(parent=self.window) or var.get()))
                browse.pack(side="right", padx=(5, 0))
                self.inputs.extend([entry, browse])
            ttk.Label(locations, text="A local source overrides the package release. Published releases use the matching "
                      "TopoPPI version when no source is found.", wraplength=650).pack(anchor="w", pady=5)
            self.log = scrolledtext.ScrolledText(output, height=8, wrap="word", state="disabled")
            self.log.pack(fill="both", expand=True)
            self.bar = ttk.Progressbar(frame, mode="indeterminate")
            self.bar.grid(row=2, column=0, sticky="ew", pady=(8, 2))
            self.status = tk.StringVar(value="Review the selected tasks, then choose the setup action below.")
            ttk.Label(frame, textvariable=self.status, wraplength=690).grid(row=3, column=0, sticky="w", pady=5)
            buttons = ttk.Frame(frame)
            buttons.grid(row=4, column=0, sticky="ew")
            buttons.columnconfigure(2, weight=1)
            self.start = ttk.Button(buttons, text=self.start_label(), command=self.run)
            self.start.grid(row=0, column=0, sticky="w")
            self.stop = ttk.Button(buttons, text="Cancel preparation", command=self.cancel_run, state="disabled")
            self.stop.grid(row=0, column=1, padx=8, sticky="w")
            self.restart = ttk.Button(buttons, text="Restart and open Affinity" if self.restart_required else "Open Affinity",
                                      command=self.continue_app, state="normal" if self.restart_required or self.dependencies_installed else "disabled")
            self.restart.grid(row=1, column=0, columnspan=2, sticky="w", pady=(6, 0))
            ttk.Button(buttons, text="Close", command=self.close).grid(row=1, column=2, sticky="e", pady=(6, 0))
            for container in (frame, choice, locations):
                container.bind("<Configure>", lambda event: [child.configure(wraplength=max(260, event.width - 24))
                    for child in event.widget.winfo_children() if child.winfo_class() == "TLabel"])
            self.prepare_models.trace_add("write", lambda *_: self.start.config(text=self.start_label()))
            for variable in self.task_vars.values():
                variable.trace_add("write", lambda *_: self.start.config(text=self.start_label()))
            self.window.protocol("WM_DELETE_WINDOW", self.close)
            if self.restart_required:
                self.status.set("Dependencies changed. Save any current figures, then restart; incomplete resource preparation can be retried first.")

        def start_label(self):
            if not self.prepare_models.get():
                return "Install dependencies"
            return "Install and prepare" if self.needs_installation() else "Prepare selected tasks"

        def needs_installation(self):
            neural_missing = any(task in {"known", "unseen"} for task in self.selected_tasks()) and any(
                importlib.util.find_spec(module) is None for module in ("torch", "esm", "transformers", "torch_geometric"))
            return not self.dependencies_installed or neural_missing or (
                self.installed_source != (self.source.get().strip() or None))

        def selected_tasks(self):
            return tuple(task for task in TASKS if self.task_vars[task].get())

        def application_busy(self):
            owner = getattr(host, "host", None) or host
            workers = (getattr(owner, "_worker_thread", None),
                       getattr(getattr(owner, "_affinity_panel", None), "worker", None), getattr(host, "worker", None))
            return getattr(owner, "_busy", False) is True or any(worker and worker.is_alive() for worker in workers)

        def continue_app(self):
            if self.restart_required:
                self.restart_app()
            elif not self.application_busy() and not (self.worker and self.worker.is_alive()):
                self.close()
                owner = getattr(host, "host", None) or host
                if hasattr(owner, "_open_affinity"):
                    owner._open_affinity()
                    panel = getattr(owner, "_affinity_panel", None)
                    if panel is not None and hasattr(panel, "restore_setup") and self.selected_tasks():
                        panel.restore_setup({"task": self.selected_tasks()[0]})

        def restart_app(self):
            if self.application_busy() or (self.worker and self.worker.is_alive()):
                self.status.set("Wait for or cancel the active task, save your work, then restart TopoPPI.")
                return
            try:
                subprocess.Popen(restart_command(host, task=next(iter(self.selected_tasks()), None)))
            except OSError as exc:
                self.status.set(f"Could not restart: {exc}. Close and reopen TopoPPI manually.")
                return
            self.close()
            owner = getattr(host, "host", None) or host
            owner.close()

        def run(self):
            if self.worker and self.worker.is_alive():
                return
            if self.application_busy():
                self.status.set("Wait for or cancel the active analysis before preparing prediction dependencies.")
                return
            download_now = self.prepare_models.get()
            chosen = self.selected_tasks()
            if download_now and not chosen:
                self.status.set("Choose at least one prediction task, or clear resource preparation to install dependencies only.")
                self.tabs.select(0)
                return
            source = self.source.get().strip() or None
            cache_dir = self.cache.get().strip() or None
            install_packages = self.needs_installation() or not download_now
            if install_packages:
                try:
                    plugin_requirement(source)
                except ValueError as exc:
                    self.status.set(str(exc))
                    self.tabs.select(1)
                    return
            self.cancel.clear()
            self.finished = False
            self.ready = False
            self.start.config(state="disabled")
            self.stop.config(state="normal")
            self.restart.config(state="disabled")
            for widget in self.inputs:
                widget.config(state="disabled")
            self.bar.start(15)
            self.tabs.select(self.progress_tab)
            self.status.set("Installing prediction dependencies..." if install_packages else "Preparing selected prediction resources...")

            def work():
                try:
                    install_affinity(source=source, download=download_now, tasks=chosen, cache_dir=cache_dir,
                                     install_packages=install_packages,
                                     progress=lambda text: self.messages.put(("log", text)), cancel_event=self.cancel,
                                     on_install_started=lambda: self.messages.put(("install_started", "Installing dependencies; a restart will apply changes.")),
                                     on_installed=lambda: self.messages.put(("installed", source)))
                except InstallationCancelled as exc:
                    self.messages.put(("cancelled", str(exc)))
                except Exception as exc:
                    self.messages.put(("error", str(exc)))
                else:
                    self.messages.put(("done", "Selected prediction resources are ready." if download_now else
                                       "Plugin installed. Prediction resources can be prepared when needed."))

            self.worker = threading.Thread(target=work, daemon=False)
            self.worker.start()
            self.poll_id = self.window.after(100, self.poll)

        def cancel_run(self):
            self.cancel.set()
            self.status.set("Cancelling preparation; waiting for the current process to stop...")
            self.stop.config(state="disabled")

        def mark_restart_required(self):
            self.restart_required = True
            owner = getattr(host, "host", None) or host
            owner._affinity_restart_required = True
            if hasattr(owner, "_refresh_action_states"):
                owner._refresh_action_states()
            if hasattr(owner, "status_var"):
                owner.status_var.set("Save current figures, then restart TopoPPI to apply affinity dependency changes.")
            self.restart.config(text="Restart and open Affinity")

        def poll(self):
            self.poll_id = None
            while not self.messages.empty():
                kind, value = self.messages.get_nowait()
                if kind == "installed":
                    self.dependencies_installed = True
                    self.installed_source = value
                    self.mark_restart_required()
                    value = "Dependencies installed. Remaining resource preparation can be retried without reinstalling packages."
                    self.status.set("Dependencies installed; preparing selected resources..." if self.prepare_models.get() else
                                    "Dependencies installed.")
                elif kind == "install_started":
                    self.mark_restart_required()
                self.log.config(state="normal")
                self.log.insert("end", value + "\n")
                self.log.see("end")
                self.log.config(state="disabled")
                if kind in {"done", "error", "cancelled"}:
                    self.finished = True
                    self.ready = kind == "done"
                    self.bar.stop()
                    if kind == "done":
                        self.status.set(value + (" Save current figures, then restart to continue." if self.restart_required else " Open Affinity to continue."))
                    elif kind == "cancelled":
                        self.status.set("Preparation cancelled. Your choices and completed downloads are retained; retry when ready.")
                    else:
                        self.status.set("Preparation failed: " + value.splitlines()[0][:180] + " See Progress for details, then retry.")
                    self.start.config(text=self.start_label() if self.ready else "Retry preparation", state="normal")
                    self.stop.config(state="disabled")
                    self.restart.config(state="normal" if self.restart_required or self.dependencies_installed else "disabled")
                    for widget in self.inputs:
                        widget.config(state="normal")
            if self.closing and (self.worker is None or not self.worker.is_alive()):
                self.close()
            elif not self.finished or (self.worker and self.worker.is_alive()):
                self.poll_id = self.window.after(100, self.poll)

        def close(self):
            self.closing = True
            if self.worker and self.worker.is_alive():
                self.cancel_run()
                if self.poll_id is None:
                    self.poll_id = self.window.after(100, self.poll)
                return False
            if self.poll_id:
                self.window.after_cancel(self.poll_id)
                self.poll_id = None
            if self.window.winfo_exists():
                self.window.destroy()
            return True

    host._affinity_installer = Installer()
    return host._affinity_installer

if __name__ == "__main__":
    raise SystemExit(main())
