"""Display optional affinity results alongside the main interface canvas."""

import csv
import math
import tkinter as tk
from copy import deepcopy
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


class _LegacyHostMethods:
    def _init_affinity_results(self):
        self._affinity_document = None
        self._affinity_snapshot = None
        self._affinity_selected_id = None
        self.affinity_frame = ttk.Frame(self.right_frame, padding=(10, 3), style="Workspace.TFrame")
        self.var_show_affinity = tk.BooleanVar(value=False)
        heading = ttk.Frame(self.affinity_frame, style="Workspace.TFrame")
        heading.pack(fill=tk.X)
        ttk.Checkbutton(heading, text="Affinity results", variable=self.var_show_affinity,
                        command=self._toggle_affinity_results).pack(side=tk.LEFT)
        self.btn_open_affinity_result = ttk.Button(heading, text="Open result...", style="Tool.TButton",
                                                  command=self.open_affinity_result)
        self.btn_open_affinity_result.pack(side=tk.RIGHT)
        self.affinity_body = ttk.Frame(self.affinity_frame, style="Workspace.TFrame")
        self.affinity_summary_var = tk.StringVar(value="Choose Binding analysis to predict, or open a saved result.")
        self.affinity_detail_var = tk.StringVar()
        summary = ttk.Label(self.affinity_body, textvariable=self.affinity_summary_var,
                            style="Workspace.TLabel", wraplength=600)
        summary.pack(fill=tk.X, pady=(3, 4))
        table_frame = ttk.Frame(self.affinity_body)
        table_frame.pack(fill=tk.X)
        self.affinity_table = ttk.Treeview(table_frame, columns=("variant", "mutations", "value", "status"),
                                           show="headings", height=3, selectmode="browse")
        for key, title, width in (("variant", "Variant", 110), ("mutations", "Substitutions", 160),
                                  ("value", "Prediction", 100), ("status", "Status", 110)):
            self.affinity_table.heading(key, text=title)
            self.affinity_table.column(key, width=width, minwidth=65)
        self.affinity_table.pack(side=tk.LEFT, fill=tk.X, expand=True)
        scroll = ttk.Scrollbar(table_frame, command=self.affinity_table.yview)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.affinity_table.configure(yscrollcommand=scroll.set)
        self.affinity_table.bind("<<TreeviewSelect>>", self._affinity_row_selected)
        controls = ttk.Frame(self.affinity_body)
        controls.pack(fill=tk.X, pady=4)
        self.affinity_partner_var = tk.StringVar()
        self.affinity_partner = ttk.Combobox(controls, textvariable=self.affinity_partner_var,
                                             state="disabled", width=13)
        self.affinity_partner.grid(row=0, column=0, sticky="ew", padx=(0, 3))
        self.affinity_partner.bind("<<ComboboxSelected>>", lambda _e: self.show_affinity_interface())
        self.btn_show_affinity = ttk.Button(controls, text="Show interface", command=self.show_affinity_interface)
        self.btn_show_affinity.grid(row=0, column=1, sticky="ew", padx=3)
        self.btn_affinity_workspace = ttk.Button(controls, text="Prediction workspace...",
                                                command=self.open_affinity_workspace)
        self.btn_affinity_workspace.grid(row=0, column=2, sticky="ew", padx=3)
        self.btn_export_affinity = ttk.Button(controls, text="Export CSV...", command=self.export_affinity_predictions)
        self.btn_export_affinity.grid(row=0, column=3, sticky="ew", padx=(3, 0))
        controls.columnconfigure(tuple(range(4)), weight=1)
        detail = ttk.Label(self.affinity_body, textvariable=self.affinity_detail_var,
                           style="Workspace.TLabel", wraplength=600)
        detail.pack(fill=tk.X, pady=(0, 4))
        self.affinity_frame.bind("<Configure>", lambda event: [
            label.configure(wraplength=max(200, event.width - 24)) for label in (summary, detail)])
        self._sync_affinity_controls()

    def _toggle_affinity_results(self):
        if self.var_show_affinity.get():
            self.affinity_body.pack(fill=tk.X)
        else:
            self.affinity_body.pack_forget()

    def open_affinity_result(self):
        if self._busy or self._closing:
            return
        path = filedialog.askopenfilename(parent=self.root, title="Open affinity analysis",
                                          filetypes=[("Affinity analysis", "*.json")])
        if not path:
            return
        try:
            from topoppi_affinity.workflow import load_analysis

            self.receive_affinity_result(load_analysis(path), show=True)
        except ImportError:
            messagebox.showinfo("Affinity results", "Install the affinity plugin through Plugins, "
                                "then open this result again.", parent=self.root)
        except (OSError, ValueError, KeyError) as exc:
            self.log(f"Cannot open affinity result: {exc}")
            messagebox.showerror("Open affinity result", str(exc), parent=self.root)

    def receive_affinity_result(self, document, *, variant_id=None, side=0, show=False):
        """Receive completed/partial results on the Tk thread without loading weights."""
        if self._closing:
            return
        report = document.report
        same = self._affinity_snapshot and self._affinity_snapshot.get("analysis_json") == str(document.path)
        if variant_id is None and same:
            variant_id = self._affinity_selected_id
            side = max(0, self.affinity_partner.current())
        snapshot = {"analysis_json": str(document.path), "analysis_id": report.get("analysis_id"),
                    "input": report.get("input", {}), "model": report.get("model", {}),
                    "groups": report.get("interface", {}).get("groups", [[], []]),
                    "predictions": document.predictions.to_dict(orient="records"),
                    "variant_id": variant_id, "side": side}
        self._affinity_document = document
        self._display_affinity_snapshot(snapshot)
        self.log(f"Affinity results loaded: {document.path}")
        if document.interface is not None and (show or self._successful_single_run is None):
            self.show_affinity_interface()

    def _display_affinity_snapshot(self, snapshot):
        """Restore numeric results embedded in an atlas without loading weights."""
        self._affinity_snapshot = deepcopy(snapshot)
        snapshot = self._affinity_snapshot
        rows = snapshot["predictions"]
        self.affinity_table.delete(*self.affinity_table.get_children())
        target = snapshot.get("model", {}).get("target")
        unit = {"pkd": "pKd", "ddg": "ddG (kcal/mol)"}.get(target, "Prediction")
        self.affinity_table.heading("value", text=unit)
        for row in rows:
            value = row.get("predicted", "")
            try:
                number = float(value)
                value = f"{number:.4g}" if math.isfinite(number) else "Unavailable"
            except (ValueError, TypeError):
                value = "Unavailable"
            status = {"ok": "Complete", "base_only": "Base model only", "error": "Failed",
                      "pending": "Pending", "cancelled": "Cancelled"}.get(row.get("status"), row.get("status", ""))
            self.affinity_table.insert("", tk.END, iid=str(row["variant_id"]),
                values=(row["variant_id"], row.get("mutations", ""), value, status))
        groups = snapshot.get("groups", [[], []])
        self.affinity_partner.configure(values=[f"{name}: {','.join(group)}" for name, group in zip(("A", "B"), groups, strict=True)])
        self.affinity_partner.current(snapshot.get("side", 0))
        available = self.affinity_table.get_children()
        chosen = snapshot.get("variant_id")
        self._affinity_selected_id = chosen if chosen in available else available[0] if available else None
        if self._affinity_selected_id:
            self.affinity_table.selection_set(self._affinity_selected_id)
            self.affinity_table.see(self._affinity_selected_id)
        snapshot["variant_id"] = self._affinity_selected_id
        source = Path(snapshot.get("input", {}).get("path") or snapshot.get("analysis_json", "Saved analysis")).name
        meaning = "Higher pKd indicates tighter binding." if target == "pkd" else (
            "Positive ddG indicates weaker binding." if target == "ddg" else "See the model details for units.")
        completed = sum(row.get("status") in {"ok", "base_only"} for row in rows)
        self.affinity_summary_var.set(f"{source} | {completed}/{len(rows)} results | {meaning}")
        self.affinity_frame.pack(fill=tk.X, before=self.canvas_frame)
        self.var_show_affinity.set(True)
        self._toggle_affinity_results()
        self._sync_affinity_controls()

    def _affinity_row_selected(self, _event=None):
        selection = self.affinity_table.selection()
        if not selection or selection[0] == self._affinity_selected_id:
            return
        self._affinity_selected_id = selection[0]
        self._affinity_snapshot["variant_id"] = selection[0]
        self._sync_affinity_controls()
        if self._affinity_document is not None and self._affinity_document.interface is not None and not self._busy:
            self.show_affinity_interface()

    def _affinity_canvas_matches(self):
        run = self._successful_single_run
        link = run["manifest"].get("affinity", {}) if run else {}
        snapshot = getattr(self, "_affinity_snapshot", None)
        return bool(snapshot and link.get("analysis_json") == snapshot["analysis_json"])

    def show_affinity_interface(self, *, mode=None, style=None):
        document = self._affinity_document
        if self._busy or self._closing or document is None:
            return False
        side = max(0, self.affinity_partner.current())
        variant = self._affinity_selected_id
        if mode is None:
            mode = self._successful_single_run["style"].get("view", "surface") if self._affinity_canvas_matches() else "surface"
        try:
            patches, viz, native_style = document.view(side, mode, variant)
            if not patches:
                raise ValueError("This partner has no retained interface. Its prediction remains in the result table.")
        except (OSError, ValueError, KeyError) as exc:
            self.affinity_detail_var.set(f"View: {exc}. Open the prediction workspace to generate 2D atlases or inspect the analysis.")
            if self._successful_single_run is not None:
                self._restore_atlas_style(self._successful_single_run["style"])
            self.log(f"Affinity view: {exc}")
            return False
        if style is None and self._affinity_canvas_matches():
            previous = self._successful_single_run
            if previous["manifest"]["affinity"]["side"] == side:
                self._capture_surface_camera()
                style = deepcopy(previous["style"])
        rendered_style = dict(style or native_style)
        rendered_style.update(view=mode, highlight_residues=native_style["highlight_residues"])
        snapshot = deepcopy(self._affinity_snapshot)
        snapshot.update(variant_id=variant, side=side)
        groups = snapshot["groups"]
        params = {"path": snapshot.get("input", {}).get("path") or str(document.path),
                  "chain_a": ",".join(groups[side]), "chain_b": ",".join(groups[1-side]),
                  "output_dir": str(document.path.parent), "auto_save": False, "min_points": 1}
        manifest = {"config": {}, "interface": document.interface.metadata(), "affinity": snapshot}
        self._pending_single_run = {"viz": viz, "patches": patches, "all_patches": patches,
                                    "params": params, "manifest": manifest, "style": rendered_style,
                                    "geometry": document.interface}
        previous = self._successful_single_run
        displayed = previous["manifest"].get("affinity", {}) if previous else {}
        same_selection = bool(displayed and all(displayed.get(key) == snapshot.get(key)
                                               for key in ("analysis_json", "side", "variant_id")))
        success = self.update_plot(viz, patches, rendered_style, complete_task=True,
                                   run_params=params, run_manifest=manifest, preserve_selection=same_selection)
        if success:
            self._sync_affinity_controls()
        return success

    def _sync_affinity_controls(self):
        if not hasattr(self, "affinity_table"):
            return
        available = self._affinity_snapshot is not None
        idle = not self._busy and not self._closing
        viewable = available and self._affinity_document is not None and self._affinity_document.interface is not None
        self.btn_open_affinity_result.configure(state=tk.NORMAL if idle else tk.DISABLED)
        self.btn_show_affinity.configure(state=tk.NORMAL if idle and viewable else tk.DISABLED)
        self.affinity_partner.configure(state="readonly" if idle and viewable else tk.DISABLED)
        self.btn_affinity_workspace.configure(state=tk.NORMAL if idle and available else tk.DISABLED)
        self.btn_export_affinity.configure(state=tk.NORMAL if available else tk.DISABLED)
        if not available:
            return
        row = next((r for r in self._affinity_snapshot["predictions"] if r["variant_id"] == self._affinity_selected_id), {})
        detail = ""
        if row.get("status") == "base_only":
            detail += "Base model only; the mutation prior was unavailable. "
        if row.get("error"):
            detail += str(row["error"]) + " "
        if self._affinity_canvas_matches():
            link = self._successful_single_run["manifest"]["affinity"]
            detail += f"Canvas: {link['variant_id']}, partner {'AB'[link['side']]}. Highlighted sites belong to this variant."
        elif self._affinity_document is None:
            detail += "Saved predictions are embedded in this atlas. Open the analysis for other variants and partners."
        elif self._affinity_document.interface is None:
            detail += "Open Prediction workspace to review the saved error and retry interface preparation."
        else:
            detail += "Choose Show interface to display the selected prediction's structure."
        run = self._successful_single_run
        if self._affinity_canvas_matches() and any(
            "uv" not in patch.metadata and "uv_global" not in patch.metadata
            for patch in run.get("all_patches", run["patches"])
        ):
            detail += " Generate 2D atlases in Prediction workspace to save an editable atlas."
        self.affinity_detail_var.set(detail)

    def open_affinity_workspace(self):
        if self._busy or self._closing or not self._affinity_snapshot:
            return
        self._open_affinity()
        panel = getattr(self, "_affinity_panel", None)
        if panel is None or not panel.window.winfo_exists():
            return
        path = Path(self._affinity_snapshot["analysis_json"])
        if path.is_file():
            panel.load(path)
            variant = self._affinity_selected_id
            if variant and variant in panel.table.get_children():
                panel.table.selection_set(variant)
            panel.side.set(str(max(0, self.affinity_partner.current())))
            panel.draw()
        else:
            panel.errors.set("Locate the saved analysis.json with Open saved result. The atlas retains its prediction table.")

    def export_affinity_predictions(self):
        if not self._affinity_snapshot:
            return
        path = filedialog.asksaveasfilename(parent=self.root, title="Export affinity predictions",
                                           initialfile="predictions.csv", defaultextension=".csv",
                                           filetypes=[("CSV", "*.csv")])
        if path:
            rows = self._affinity_snapshot["predictions"]
            try:
                with open(path, "w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
                    writer.writeheader()
                    writer.writerows(rows)
                self.log(f"Affinity predictions exported to {path}")
            except OSError as exc:
                messagebox.showerror("Export predictions", str(exc), parent=self.root)





    def _open_affinity(self):
        if self._busy or self._closing:
            return
        if getattr(self, "_affinity_restart_required", False):
            from topoppi_affinity.installer import open_installer
            return open_installer(self)
        try:
            from topoppi_affinity.gui import open_panel
            return open_panel(self)
        except ImportError:
            from topoppi_affinity.installer import open_installer
            return open_installer(self)


class AffinityController:
    """Task-owned UI using the Core's small GUI extension lifecycle."""

    def __init__(self, host):
        self.host = host
        self.menu_entries = []

    def initialize(self):
        from types import MethodType

        # Explicit legacy host API for existing callers; Core knows no task names.
        for name, method in vars(_LegacyHostMethods).items():
            if callable(method):
                setattr(self.host, name, MethodType(method, self.host))
        host = self.host
        host._init_affinity_results()
        if hasattr(host, "register_binding_task"):
            for task, title, description in (
                ("absolute", "Binding affinity", "Predict whole-complex pKd and Kd."),
                ("mutation", "Mutation effects", "Predict ddG for individual or joint substitutions."),
            ):
                host.register_binding_task(task, title=title, description=description, plugin="affinity",
                    open_task=lambda setup=None, task=task: self.open(setup={**(setup or {}), "prediction_task": task}))
            host.analysis_menu.add_command(label="Open affinity result...", command=host.open_affinity_result)
            self.menu_entries.append(host.analysis_menu.index("end"))
            self._shared_navigation = True
        else:
            for label, callback in (("Predict affinity...", host._open_affinity),
                                    ("Open affinity result...", host.open_affinity_result)):
                host.analysis_menu.add_command(label=label, command=callback)
                self.menu_entries.append(host.analysis_menu.index("end"))
            host.btn_affinity = ttk.Button(host.plugin_actions_frame, text="Affinity", width=8, style="Tool.TButton",
                                           command=host._open_affinity)
            host.btn_affinity.pack(side=tk.RIGHT)
            self._shared_navigation = False

    def open(self, setup=None):
        existing = getattr(self.host, "_affinity_panel", None)
        had_existing = bool((existing is not None and existing.window.winfo_exists())
                            or getattr(self.host, "_affinity_draft", None))
        if setup is not None and existing is not None and existing.window.winfo_exists() and (
                existing._working or existing.closing or (existing.worker is not None and existing.worker.is_alive())):
            raise RuntimeError("A binding prediction is active. Wait for it to finish or cancel it before changing tasks or adding variants.")
        panel = self.host._open_affinity()
        if setup is not None:
            panel = panel or getattr(self.host, "_affinity_panel", None)
            if panel is None or not hasattr(panel, "restore_setup"):
                raise RuntimeError("The affinity prediction editor is not ready. Complete dependency setup and restart TopoPPI before opening this task.")
            if not had_existing and any(key in setup for key in ("structure", "chains_a", "chains_b")):
                panel.restore_setup(setup, initial=True)
            else:
                panel.restore_setup(setup)
        return panel

    @property
    def restart_required(self):
        return bool(getattr(self.host, "_affinity_restart_required", False))

    def sync_controls(self):
        host = self.host
        state = tk.DISABLED if host._busy or host._closing or self.restart_required else tk.NORMAL
        if hasattr(host, "btn_affinity"):
            host.btn_affinity.configure(state=state)
        for position, index in enumerate(self.menu_entries):
            menu_state = state if position == 0 and not getattr(self, "_shared_navigation", False) else tk.DISABLED if host._busy or host._closing else tk.NORMAL
            host.analysis_menu.entryconfigure(index, state=menu_state)
        host._sync_affinity_controls()

    def redraw(self, style):
        host = self.host
        run = host._successful_single_run
        if (host._affinity_canvas_matches() and host._affinity_document is not None
                and style.get("view") != run["style"].get("view")):
            host.show_affinity_interface(mode=style.get("view"), style=style)
            return True
        return False

    def document_loaded(self, manifest):
        if snapshot := manifest.get("affinity"):
            run = self.host._successful_single_run
            if run is not None:
                run["manifest"].setdefault("interface", {
                    "source": snapshot.get("input", {}), "groups": snapshot.get("groups", []),
                })
            self.host._affinity_document = None
            self.host._display_affinity_snapshot(snapshot)

    def filename(self, run):
        import re

        if link := run["manifest"].get("affinity"):
            source = Path(str(run["params"]["path"]))
            variant = re.sub(r"[^\w.-]+", "_", str(link.get("variant_id") or "complex"))
            return f"{source.stem}_affinity_{variant}_partner_{'AB'[link['side']]}.png"

    def diagnostics(self):
        snapshot = getattr(self.host, "_affinity_snapshot", None)
        return {"analysis": snapshot.get("analysis_json") if snapshot else None}

    def close(self):
        for name in ("_affinity_installer", "_affinity_panel"):
            child = getattr(self.host, name, None)
            if child is not None and not child.close():
                return False
        return True


def create_controller(host):
    return AffinityController(host)
