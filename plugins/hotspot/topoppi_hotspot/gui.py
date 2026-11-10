"""Hotspot desktop workspace, loaded only when its GUI contribution is selected."""
from __future__ import annotations

import gc
import json
import queue
import threading
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, ttk

from .gui_workflow import RankingDocument, RankingSettings, load_ranking, run_ranking


class HotspotPanel:
    def __init__(self, controller):
        self.controller, self.host = controller, controller.host
        self.window = tk.Toplevel(self.host.root)
        self.window.title("TopoPPI Hotspot | Rank interface residues")
        self.window.geometry("1040x760")
        self.window.minsize(850, 650)
        self.window.protocol("WM_DELETE_WINDOW", self.close)
        self.vars, self.inputs = {}, []
        self.document = None
        self.worker = None
        self._working = False
        self.closing = False
        self.cancel = threading.Event()
        self.events = queue.Queue()
        self.status = tk.StringVar(value="Choose a partner and results folder to predict alanine-substitution binding changes.")
        self.error = tk.StringVar()
        self.model_status = tk.StringVar(value="Mutation-derived surface model · CPU · requires residue-aware OptCuts.")
        self.summary = tk.StringVar(value="Predicts Ala ddG at standard non-alanine interface residues. Larger positive values indicate greater loss of binding.")
        self.selection_status = tk.StringVar(value="Select ranked residues to create variants.")
        outer = ttk.Frame(self.window, padding=10)
        outer.pack(fill=tk.BOTH, expand=True)
        navigation = getattr(self.host, "add_binding_navigation", None)
        if navigation is not None:
            navigation(outer, current_task="hotspot", setup_provider=self.binding_setup)
        header = ttk.Frame(outer)
        header.pack(fill=tk.X, pady=(0, 4))
        ttk.Label(header, text="Rank interface residues", font=("TkDefaultFont", 15, "bold")).pack(side=tk.LEFT)
        self.open_button = ttk.Button(header, text="Open saved ranking...", command=self.open_saved)
        self.open_button.pack(side=tk.RIGHT)

        setup = ttk.LabelFrame(outer, text="Structure and model", padding=6)
        setup.pack(fill=tk.X)
        setup.columnconfigure(1, weight=1)
        self._path_row(setup, 0, "structure", "Structure file", "file")
        groups = ttk.Frame(setup)
        groups.grid(row=1, column=0, columnspan=3, sticky="ew", pady=3)
        for column, (key, label) in enumerate((("chains_a", "Partner 1 chains"), ("chains_b", "Partner 2 chains"))):
            ttk.Label(groups, text=label).grid(row=0, column=column * 2, sticky="w", padx=(0, 6))
            self.vars[key] = tk.StringVar()
            entry = ttk.Entry(groups, textvariable=self.vars[key], width=12)
            entry.grid(row=0, column=column * 2 + 1, sticky="ew", padx=(0, 18))
            self.inputs.append((entry, "normal"))
        ttk.Label(groups, text="Rank").grid(row=0, column=4, sticky="w", padx=(0, 6))
        self.partner = tk.StringVar(value="Partner 1")
        choice = ttk.Combobox(groups, textvariable=self.partner, values=("Partner 1", "Partner 2"), state="readonly", width=12)
        choice.grid(row=0, column=5, sticky="w")
        self.inputs.append((choice, "readonly"))
        ttk.Label(setup, text="Use author chain IDs. Separate chains in one partner with commas, for example A,B.").grid(
            row=2, column=1, columnspan=2, sticky="w", pady=(0, 2))
        self._path_row(setup, 3, "model_dir", "Custom model (optional)", "folder")
        ttk.Label(setup, textvariable=self.model_status, wraplength=650).grid(row=4, column=1, columnspan=2, sticky="w", pady=(0, 2))
        self._path_row(setup, 5, "output_dir", "Results folder", "folder")
        controls = ttk.Frame(setup)
        controls.grid(row=6, column=1, columnspan=2, sticky="ew", pady=(4, 0))
        for title, command in (("Use main selection", self.use_main_selection), ("New results folder", self.new_output)):
            button = ttk.Button(controls, text=title, command=command)
            button.pack(side=tk.LEFT, padx=(0, 6))
            self.inputs.append((button, "normal"))
        self.sequence_window = tk.Toplevel(self.window)
        self.sequence_window.withdraw()
        self.sequence_window.title("Hotspot · Sequence model options")
        self.sequence_window.geometry("760x160")
        self.sequence_window.minsize(700, 150)
        self.sequence_window.transient(self.window)
        self.sequence_window.protocol("WM_DELETE_WINDOW", self.sequence_window.withdraw)
        advanced = ttk.LabelFrame(self.sequence_window, text="Prepared sequence embeddings", padding=10)
        advanced.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.sequence_button = ttk.Button(controls, text="Sequence model options...", command=self.show_sequence_options)
        self.sequence_button.pack(side=tk.LEFT, padx=(0, 6))
        self.inputs.append((self.sequence_button, "normal"))
        advanced.columnconfigure(1, weight=1)
        self._path_row(advanced, 0, "embedding_registry", "Embedding registry", "folder")
        ttk.Label(advanced, text="Select prepared sequence embeddings when required by the model. Ranking runs on CPU.", wraplength=500).grid(
            row=1, column=1, columnspan=2, sticky="w")
        ttk.Button(advanced, text="Done", command=self.sequence_window.withdraw).grid(row=2, column=2, sticky="e", pady=(6, 0))

        actions = ttk.Frame(outer)
        actions.pack(fill=tk.X, pady=(0, 6))
        self.run_button = ttk.Button(actions, text="Rank residues", command=self.start)
        self.run_button.pack(side=tk.LEFT)
        self.cancel_button = ttk.Button(actions, text="Cancel", command=self.cancel_run, state=tk.DISABLED)
        self.cancel_button.pack(side=tk.LEFT, padx=6)
        self.progress = ttk.Progressbar(actions, mode="indeterminate")
        status_label = ttk.Label(outer, textvariable=self.status, wraplength=970)
        status_label.pack(fill=tk.X, pady=2)
        error_label = ttk.Label(outer, textvariable=self.error, foreground="#a32f2f", wraplength=970)
        outer.bind("<Configure>", lambda event: [label.configure(wraplength=max(1, event.width - 24))
                                                  for label in (status_label, error_label)])
        results = ttk.LabelFrame(outer, text="Ranked residues", padding=8)
        results.pack(fill=tk.BOTH, expand=True)
        self.error.trace_add("write", lambda *_: error_label.pack(fill=tk.X, pady=(0, 6), before=results)
                             if self.error.get() else error_label.pack_forget())
        results.columnconfigure(0, weight=1)
        results.rowconfigure(1, weight=1)
        summary_label = ttk.Label(results, textvariable=self.summary, wraplength=780)
        summary_label.grid(row=0, column=0, sticky="ew", pady=(0, 4))
        table_frame = ttk.Frame(results)
        table_frame.grid(row=1, column=0, sticky="nsew")
        self.table = ttk.Treeview(table_frame, columns=("rank", "residue", "score", "percentile", "flags"), show="headings", selectmode="extended", height=6)
        for key, label, width in (("rank", "Rank", 65), ("residue", "Residue", 140), ("score", "ddG (kcal/mol)", 115),
                                  ("percentile", "Percentile", 100), ("flags", "Coverage", 270)):
            self.table.heading(key, text=label)
            self.table.column(key, width=width, minwidth=50)
        self.table.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.table.yview)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.table.configure(yscrollcommand=scroll.set)
        self.table.bind("<<TreeviewSelect>>", self.selected)
        selection_label = ttk.Label(results, textvariable=self.selection_status, wraplength=780)
        selection_label.grid(row=2, column=0, sticky="ew", pady=(4, 0))
        review = ttk.Frame(results)
        review.grid(row=3, column=0, sticky="ew", pady=(4, 0))
        results.bind("<Configure>", lambda event: [label.configure(wraplength=max(1, event.width - 24))
                                                    for label in (summary_label, selection_label)])
        self.variants_button = ttk.Button(review, text="Create variants...", command=self.create_variants, state=tk.DISABLED)
        self.variants_button.pack(side=tk.LEFT, padx=(0, 6))
        self.review_buttons = []
        for title, callback in (("Show ranked surface", self.show_surface), ("Export rankings...", self.export),
                                ("Export residue values...", lambda: self.export(annotations=True))):
            button = ttk.Button(review, text=title, command=callback, state=tk.DISABLED)
            button.pack(side=tk.LEFT, padx=(0, 6))
            self.review_buttons.append(button)
        self.use_main_selection()
        self.vars["model_dir"].trace_add("write", lambda *_: self.inspect_model())
        self._poll_id = self.window.after(80, self._poll)

    def show_sequence_options(self):
        self.sequence_window.deiconify()
        self.sequence_window.lift()

    def _path_row(self, parent, row, key, label, kind):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 12), pady=2)
        self.vars[key] = tk.StringVar()
        entry = ttk.Entry(parent, textvariable=self.vars[key])
        entry.grid(row=row, column=1, sticky="ew", pady=2)
        button = ttk.Button(parent, text="Browse...", command=lambda: self.browse(key, kind))
        button.grid(row=row, column=2, padx=(8, 0), pady=2)
        self.inputs.extend(((entry, "normal"), (button, "normal")))

    def browse(self, key, kind):
        value = (filedialog.askdirectory(parent=self.window, title="Choose " + key.replace("_", " ")) if kind == "folder"
                 else filedialog.askopenfilename(parent=self.window, title="Choose structure", filetypes=[("Structures", "*.pdb *.cif *.mmcif"), ("All files", "*")]))
        if value:
            self.vars[key].set(value)

    def use_main_selection(self):
        for key, name in (("structure", "var_input_path"), ("chains_a", "var_chain_a"), ("chains_b", "var_chain_b")):
            self.vars[key].set(getattr(self.host, name).get())
        if not self.vars["output_dir"].get():
            self.new_output()

    def new_output(self):
        source = Path(self.vars["structure"].get()).expanduser()
        parent = source.parent if source.is_file() else Path.home()
        self.vars["output_dir"].set(str(parent / "topoppi-results" / ("hotspot-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))))

    def inspect_model(self):
        from .resources import resolve_model_dir

        try:
            path = resolve_model_dir(self.vars["model_dir"].get()) / "model.json"
            metadata = json.loads(path.read_text(encoding="utf-8"))
            label = metadata.get("representation") or path.parent.name
            self.model_status.set(f"{label} · {metadata['kind']} model. Uses this model's feature settings.")
        except (OSError, ValueError, KeyError, TypeError):
            self.model_status.set("Choose a folder containing model.json, or clear the field to use the official model.")

    def start(self):
        if self._working or self.host._busy or self.host._closing:
            return
        values = {key: variable.get().strip() for key, variable in self.vars.items()}
        from .resources import resolve_model_dir

        try:
            resolve_model_dir(values["model_dir"])
        except ValueError as exc:
            self.error.set(str(exc))
            return
        if not values["output_dir"]:
            self.error.set("Choose a results folder, or select New results folder.")
            return
        partner = values["chains_a"] if self.partner.get() == "Partner 1" else values["chains_b"]
        settings = RankingSettings(**values, partner=partner)
        getter = getattr(self.host, "get_current_geometry", None)
        prepared = getter() if getter else None
        self.error.set("")
        self.cancel.clear()
        self._set_working(True)
        self.status.set("Preparing the ranking with the selected model...")

        def work():
            from topoppi import ExecutionContext
            from topoppi.errors import ExecutionCancelled

            try:
                context = ExecutionContext(cancel_event=self.cancel, progress=lambda event: self.events.put(("progress", event)))
                result = run_ranking(settings, prepared=prepared, context=context)
                self.events.put(("result", result))
            except (ExecutionCancelled, InterruptedError):
                self.events.put(("cancelled", "Ranking cancelled. Choose a new results folder to run again."))
            except Exception as exc:
                self.events.put(("error", str(exc)))
            finally:
                self.events.put(("done", None))

        # Collect discarded host plots on the UI thread before computation can
        # trigger their Tk and PIL image finalizers in the ranking worker.
        gc.collect()
        self.worker = threading.Thread(target=work, name="topoppi-hotspot", daemon=False)
        self.worker.start()

    def _set_working(self, value):
        self._working = value
        for widget, state in self.inputs:
            widget.configure(state=tk.DISABLED if value else state)
        self.run_button.configure(state=tk.DISABLED if value else tk.NORMAL)
        self.open_button.configure(state=tk.DISABLED if value else tk.NORMAL)
        self.cancel_button.configure(state=tk.NORMAL if value else tk.DISABLED)
        if value:
            self.progress.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=8)
            self.progress.start()
        else:
            self.progress.stop()
            self.progress.pack_forget()
        for button in self.review_buttons:
            button.configure(state=tk.NORMAL if self.document is not None and not value else tk.DISABLED)
        self.selected()

    def cancel_run(self):
        self.cancel.set()
        self.status.set("Stopping after the current calculation finishes...")
        self.cancel_button.configure(state=tk.DISABLED)

    def _poll(self):
        self._poll_id = None
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "progress" and value.get("message"):
                    self.status.set(value["message"])
                elif kind == "result" and not self.closing:
                    self.display(value)
                elif kind == "error":
                    self.error.set(value)
                    self.status.set("Review the inputs, then run the ranking again.")
                elif kind == "cancelled":
                    self.status.set(value)
                elif kind == "done":
                    self._set_working(False)
        except queue.Empty:
            pass
        if self.closing and self.close():
            return
        self._poll_id = self.window.after(80, self._poll)

    def display(self, document, selected=()):
        self.document = document
        self.controller.document = document
        self.table.delete(*self.table.get_children())
        labels = {"missing_native_footprint": "Native footprint unavailable", "missing_CA": "CA unavailable",
                  "missing_spatial_context": "Spatial context imputed", "tied_score": "Tied score"}
        is_ddg = document.metadata.get("target") == "ddg"
        self.table.heading("score", text="ddG (kcal/mol)" if is_ddg else "Score")
        for row in document.rows:
            low, high = int(row["rank_min"]), int(row["rank_max"])
            rank = str(low) if low == high else f"{low}-{high}"
            percentile = f"{float(row.get('full_candidate_percentile', 0)):.1%}"
            flags = "; ".join(labels.get(flag, flag) for flag in str(row.get("flags", "")).split(";") if flag) or "Complete"
            self.table.insert("", tk.END, iid=row["residue_key"], values=(rank, row["residue_key"], f"{float(row['rank_score']):.5g}", percentile, flags))
        if document.rows:
            selected_keys = [key for key in selected if self.table.exists(key)]
            self.table.selection_set(selected_keys or [document.rows[0]["residue_key"]])
        partner = ",".join(document.metadata.get("selected_partner", []))
        source = Path(document.metadata.get("structure", "Saved structure")).name
        temperature = float(document.metadata.get("temperature_K", 298.15))
        meaning = (f"Predicted Ala ddG (kcal/mol) at {temperature:g} K. Larger positive values indicate greater loss of binding."
                   if is_ddg else "Standard non-alanine interface residues. Higher scores and within-partner percentiles rank first.")
        self.summary.set(f"{source} · partner {partner} · {len(document.rows)} ranked residues.\n{meaning}")
        self.status.set(f"Ranking ready · {document.directory}")
        self.error.set("")
        for button in self.review_buttons:
            button.configure(state=tk.NORMAL)
        self.selected()

    def selected(self, _event=None):
        selection = self.table.selection()
        self.controller.selected = selection[0] if selection else None
        self.controller.selected_residues = selection
        count = len(selection)
        self.selection_status.set(
            f"{count} residue{'s' if count != 1 else ''} selected. Ctrl/Cmd-click or Shift-click to select more. "
            "Variant creation uses the displayed ranking's structure and partners."
            if count else "Select ranked residues to create variants.")
        self.variants_button.configure(state=tk.NORMAL if count and self.document is not None and not self._working else tk.DISABLED)

    def binding_setup(self, displayed=False):
        """Navigate with current inputs or explicitly continue the displayed result."""
        if not displayed or self.document is None:
            return {name: self.vars[name].get() for name in ("structure", "chains_a", "chains_b")}
        metadata = self.document.metadata
        groups = metadata.get("partner_groups", [])
        if not metadata.get("structure") or len(groups) != 2:
            raise ValueError("This ranking does not record its source structure and both partner groups. Open a complete ranking folder.")
        return {"structure": metadata["structure"], "chains_a": ",".join(groups[0]), "chains_b": ",".join(groups[1])}

    def create_variants(self):
        if self.document is None or not self.table.selection():
            self.error.set("Select one or more ranked residues to create variants.")
            return False
        if self._working or self.host._busy or self.host._closing:
            self.error.set("Wait for the current calculation to finish before creating variants.")
            return False
        try:
            create = getattr(self.host, "create_variants_from_residues", None)
            if create is None:
                raise ValueError("Open Hotspot from the Binding analysis workspace to create variants with Affinity.")
            selection = set(self.table.selection())
            residues = []
            for row in self.document.rows:
                if row["residue_key"] not in selection:
                    continue
                if any(name not in row for name in ("chain", "residue_number", "insertion_code", "wt")):
                    raise ValueError("This saved ranking lacks explicit residue identities. Rerun the ranking to create variants; its table can still be reviewed and exported.")
                residues.append({"chain": row["chain"], "position": int(row["residue_number"]),
                                 "insertion": row["insertion_code"], "wt": row["wt"], "residue_key": row["residue_key"]})
            result = create(residues, setup=self.binding_setup(displayed=True))
            self.error.set("")
            return result
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            self.error.set(str(exc))
            return False

    def open_saved(self):
        path = filedialog.askopenfilename(parent=self.window, title="Open Hotspot ranking", filetypes=[("Ranking metadata", "metadata.json"), ("JSON", "*.json")])
        if path:
            try:
                self.display(load_ranking(path))
            except (OSError, ValueError, KeyError) as exc:
                self.error.set(f"Open the folder's metadata.json with its rankings.csv. Details: {exc}")

    def show_surface(self):
        if self.document is None:
            return False
        try:
            return self.controller.show_surface(self.document, self.table.selection())
        except (OSError, ValueError, KeyError) as exc:
            self.error.set(str(exc))
            return False

    def export(self, *, annotations=False):
        if self.document is None:
            return
        path = filedialog.asksaveasfilename(parent=self.window, title="Export residue scores" if annotations else "Export rankings",
            initialfile="hotspot-scores.csv" if annotations else "rankings.csv", defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if path:
            try:
                self.document.export(path, annotations=annotations)
                self.status.set(f"CSV saved · {path}")
                self.error.set("")
            except OSError as exc:
                self.error.set(f"Choose a writable CSV destination. Details: {exc}")

    def close(self):
        if not self.window.winfo_exists():
            return True
        self.closing = True
        self.cancel.set()
        if self.worker is not None and self.worker.is_alive():
            self.status.set("Closing after the current calculation finishes...")
            return False
        self.controller.draft = {"fields": {name: variable.get() for name, variable in self.vars.items()},
                                 "partner": self.partner.get()}
        if self._poll_id is not None:
            self.window.after_cancel(self._poll_id)
            self._poll_id = None
        self.window.destroy()
        return True


class HotspotController:
    def __init__(self, host):
        self.host, self.panel, self.document, self.selected = host, None, None, None
        self.selected_residues = ()
        self.draft = None
        self.button, self.menu_index = None, None

    def initialize(self):
        register = getattr(self.host, "register_binding_task", None)
        if register is not None:
            register("hotspot", title="Rank hotspots", description="Rank interface residues within one binding partner.",
                     plugin="hotspot", open_task=self.open)
            return
        self.host.analysis_menu.add_command(label="Rank hotspots...", command=self.open)
        self.menu_index = self.host.analysis_menu.index("end")
        self.button = ttk.Button(self.host.plugin_actions_frame, text="Hotspot", width=8, style="Tool.TButton", command=self.open)
        self.button.pack(side=tk.RIGHT, padx=4)

    def open(self, setup=None):
        if self.host._busy or self.host._closing:
            return None
        if self.panel is None or not self.panel.window.winfo_exists():
            self.panel = HotspotPanel(self)
            if self.draft is not None:
                for name, value in self.draft["fields"].items():
                    self.panel.vars[name].set(value)
                self.panel.partner.set(self.draft["partner"])
            if self.document is not None:
                self.panel.display(self.document, selected=self.selected_residues)
            if self.draft is not None:
                self.panel.status.set("Ranking setup restored. Review the inputs before running; the displayed ranking is unchanged."
                                      if self.document is not None else "Ranking setup restored. Review the inputs, then select Rank residues.")
        else:
            self.panel.window.deiconify()
            self.panel.window.lift()
        if setup:
            if self.panel._working:
                raise ValueError("Wait for the Hotspot ranking to finish before changing its inputs.")
            for name, value in setup.items():
                if name in self.panel.vars:
                    self.panel.vars[name].set(str(value))
        return self.panel

    def sync_controls(self):
        state = tk.DISABLED if self.host._busy or self.host._closing else tk.NORMAL
        if self.button is not None:
            self.button.configure(state=state)
        if self.menu_index is not None:
            self.host.analysis_menu.entryconfigure(self.menu_index, state=state)

    def show_surface(self, document, selected=None):
        if self.host._busy or self.host._closing:
            return False
        patches, visualizer, style = document.view(selected)
        selected_keys = [selected] if isinstance(selected, str) else list(selected or ())
        group = document.metadata["selected_partner"]
        other = next(item for item in document.metadata["partner_groups"] if item != group)
        params = {"path": document.metadata["structure"], "chain_a": ",".join(group), "chain_b": ",".join(other),
                  "output_dir": str(document.directory), "auto_save": False, "min_points": 1}
        manifest = {"config": {}, "interface": document.geometry.metadata(),
                    "hotspot": {"directory": str(document.directory), "metadata": document.metadata,
                                                "rows": document.rows, "selected": selected_keys[0] if selected_keys else None,
                                                "selected_residues": selected_keys}}
        self.host._pending_single_run = {"viz": visualizer, "patches": patches, "all_patches": patches,
            "params": params, "manifest": manifest, "style": style, "geometry": document.geometry}
        return self.host.update_plot(visualizer, patches, style, complete_task=True, run_params=params, run_manifest=manifest)

    def document_loaded(self, manifest):
        if snapshot := manifest.get("hotspot"):
            metadata = snapshot["metadata"]
            if metadata.get("structure") and len(metadata.get("partner_groups", [])) == 2:
                interface = manifest.setdefault("interface", {
                    "source": {"path": metadata["structure"]}, "groups": metadata["partner_groups"],
                })
                if run := getattr(self.host, "_successful_single_run", None):
                    run["manifest"].setdefault("interface", interface)
            self.document = RankingDocument(Path(snapshot["directory"]), snapshot["metadata"], snapshot["rows"])
            self.selected = snapshot.get("selected")
            self.selected_residues = tuple(snapshot.get("selected_residues", [self.selected] if self.selected else []))
            if self.panel is not None and self.panel.window.winfo_exists():
                self.panel.display(self.document, selected=self.selected_residues)

    def filename(self, run):
        if snapshot := run["manifest"].get("hotspot"):
            return f"{Path(snapshot['metadata']['structure']).stem}_hotspot.png"

    def diagnostics(self):
        return {"result": str(self.document.directory) if self.document else None,
                "running": bool(self.panel and self.panel._working), "selected_residue": self.selected,
                "selected_residues": list(self.selected_residues),
                "inputs": {name: value.get() for name, value in self.panel.vars.items()}
                if self.panel and self.panel.window.winfo_exists() else {}}

    def close(self):
        return self.panel.close() if self.panel is not None else True


def create_controller(host):
    return HotspotController(host)
