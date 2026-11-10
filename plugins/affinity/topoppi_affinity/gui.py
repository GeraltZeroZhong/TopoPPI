"""Optional task-oriented Tk panel. Saved-result review never loads model weights."""
from __future__ import annotations

import gc
import importlib.util
import json
import math
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
import traceback
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk
from tkinter import font as tkfont

TASKS = {
    "Binding affinity (pKd)": "absolute",
    "Mutation effects (ddG)": "mutation",
}
TASK_HELP = {
    "absolute": "Predict binding affinity for the supplied complex: pKd = -log10(Kd/M).",
    "mutation": "Predict the effect of specified substitutions: ddG = G_mut - G_wt (kcal/mol). Positive values indicate weaker binding.",
}
MODEL_SOURCES = {"Prepared task model": "prepared", "Local model file": "local"}
MUTATION_MODELS = {"Bundled physical model": "", "Known-system neural ensemble": "known",
                   "Unseen-system neural ensemble": "unseen"}
SCOPE_HELP = {
    "": "The bundled physical model supports substitutions on your supplied WT complex. Advanced known/unseen ensembles use their documented training scope.",
    "known": "Known: use only when this binding system is represented in the model training data. This is not determined by having an experimental structure.",
    "unseen": "Unseen: transfer to a binding system outside the model training data. Check the model's training membership before choosing.",
}
VIEWS = {"3D interface": "surface", "2D atlas": "atlas"}
DESCRIPTORS = {"Uniform colour": "none", "Contact count": "contact_degree",
               "Footprint area": "footprint_area_A2", "Surface area": "surface_area_A2"}
AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")


def variant_lines(text):
    """Allow plain mutation rows or ID<TAB>mutation rows, preserving sample IDs."""
    import pandas as pd

    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if "\t" in line:
            identifier, mutations = line.split("\t", 1)
        else:
            identifier, mutations = f"variant_{len(rows) + 1:03d}", line
        rows.append({"variant_id": identifier.strip(), "mutations": mutations.strip()})
    if not rows:
        return None
    from topoppi_affinity.workflow import variants_table

    return variants_table(variants=pd.DataFrame(rows))[["variant_id", "mutations"]]


class AffinityPanel:
    def __init__(self, root, *, host=None, analysis=None):
        self.window = tk.Toplevel(root)
        self.root = self.window
        self.window.title("TopoPPI Affinity — predict and review")
        width = min(1240, max(800, self.window.winfo_screenwidth() - 80))
        height = min(850, max(540, self.window.winfo_screenheight() - 100))
        self.window.geometry(f"{width}x{height}")
        # The shared navigation needs one additional row above the result canvas.
        self.window.minsize(760, 540 if host is not None and hasattr(host, "add_binding_navigation") else 500)
        self.window.protocol("WM_DELETE_WINDOW", self.close)
        self.host, self.document, self.canvas, self.figure = host, None, None, None
        self.worker = None
        self._working = False
        self._input_states = {}
        self.events = queue.Queue()
        self.cancel = threading.Event()
        self.closing = False
        self.vars = {}
        self.structure_records = None
        self._drawn_selection = ()
        self.log_history = []
        self.status = tk.StringVar(value="Choose a structure to begin, or open a saved result.")
        self.summary = tk.StringVar(value="Predictions and the selected interface will appear here.")
        self.details, self.errors = tk.StringVar(), tk.StringVar()
        self._export_previous_error = ""
        self.model_status, self.resource_details, self.task_help = tk.StringVar(), tk.StringVar(), tk.StringVar()
        self.chain_status = tk.StringVar(value="Choose a PDB or mmCIF file, then assign its chains to two binding partners.")
        self.main_input_context = tk.StringVar()
        self.task = tk.StringVar(value=next(iter(TASKS)))
        self.mutation_context = tk.StringVar(value="")
        self.mutation_model_choice = tk.StringVar(value="Bundled physical model")
        self.model_source = tk.StringVar(value="Prepared task model")
        self.scope_help, self.local_model_details = tk.StringVar(), tk.StringVar()
        self._local_model_info = None
        self._local_model_path = None
        self.side, self.partner_label = tk.StringVar(value="0"), tk.StringVar(value="A")
        self.mode, self.descriptor = tk.StringVar(value="surface"), tk.StringVar(value="none")
        self.view_label, self.descriptor_label = tk.StringVar(value="3D interface"), tk.StringVar(value="Uniform colour")
        self.partner_labels = ["A", "B"]
        self.review_buttons = []
        self.resource_state = None
        if host is not None and hasattr(host, "add_binding_navigation"):
            host.add_binding_navigation(self.window, current_task=lambda: TASKS[self.task.get()], setup_provider=self.binding_setup)
        split = ttk.Panedwindow(self.window, orient=tk.HORIZONTAL)
        split.pack(fill=tk.BOTH, expand=True)
        left, right = ttk.Frame(split, padding=8), ttk.Frame(split, padding=8)
        split.add(left, weight=0)
        split.add(right, weight=1)
        setup_files = ttk.Frame(left)
        setup_files.pack(fill=tk.X, pady=(0, 5))
        self.open_setup_button = ttk.Button(setup_files, text="Open setup...", command=self.open_setup)
        self.open_setup_button.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.save_setup_button = ttk.Button(setup_files, text="Save setup...", command=self.save_setup)
        self.save_setup_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))
        self.tabs = ttk.Notebook(left)
        self.tabs.pack(fill=tk.BOTH, expand=True)
        inputs, advanced = self._scroll_tab(self.tabs, "Set up prediction"), self._scroll_tab(self.tabs, "Advanced")

        def selected(name, default=""):
            return getattr(host, name).get() if host and hasattr(host, name) else default

        structure = ttk.LabelFrame(inputs, text="1. Structure and binding partners", padding=6)
        structure.pack(fill=tk.X, pady=(0, 8))
        self._entry(structure, "structure", "Structure file (PDB / mmCIF)", selected("var_input_path"), browse="file")
        row = ttk.Frame(structure)
        row.pack(fill=tk.X, pady=4)
        ttk.Button(row, text="Read chains", command=self.read_structure).pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.example_button = ttk.Button(row, text="Try 1BVK example...", command=self.load_example)
        self.example_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))
        self.use_main_input_button = ttk.Button(structure, text="Use main window input", command=self.use_main_window_input)
        if host is not None:
            self.use_main_input_button.pack(fill=tk.X, pady=(2, 3))
            ttk.Label(structure, textvariable=self.main_input_context, wraplength=325).pack(fill=tk.X, pady=(0, 3))
        ttk.Label(structure, textvariable=self.chain_status, wraplength=325).pack(fill=tk.X, pady=3)
        self.chain_choices = tk.Listbox(structure, height=3, selectmode=tk.EXTENDED, exportselection=False)
        self.chain_choices.pack(fill=tk.X)
        chain_buttons = ttk.Frame(structure)
        chain_buttons.pack(fill=tk.X, pady=3)
        for side in ("a", "b"):
            ttk.Button(chain_buttons, text=f"Assign to partner {side.upper()}",
                       command=lambda side=side: self._assign_chains(side)).pack(side=tk.LEFT, expand=True, fill=tk.X)
        self._entry(structure, "chains_a", "Partner A chains", selected("var_chain_a", "A"))
        self._entry(structure, "chains_b", "Partner B chains", selected("var_chain_b", "B"))
        ttk.Label(structure, text="Use commas for a partner with several chains, such as A,B.", wraplength=325).pack(fill=tk.X, pady=3)

        task_box = ttk.LabelFrame(inputs, text="2. Prediction task", padding=6)
        task_box.pack(fill=tk.X)
        self.task_selector = ttk.Combobox(task_box, textvariable=self.task, values=tuple(TASKS), state="readonly", width=36)
        self.task_selector.pack(fill=tk.X)
        ttk.Label(task_box, textvariable=self.task_help, wraplength=325).pack(fill=tk.X, pady=5)
        self.scope_box = ttk.LabelFrame(advanced, text="Mutation model", padding=6)
        self.scope_box.pack(fill=tk.X, pady=(0, 8))
        self.scope_selector = ttk.Combobox(self.scope_box, textvariable=self.mutation_model_choice,
                                          values=tuple(MUTATION_MODELS), state="readonly")
        self.scope_selector.pack(fill=tk.X)
        self.scope_selector.bind("<<ComboboxSelected>>", lambda _event:
            self.mutation_context.set(MUTATION_MODELS[self.mutation_model_choice.get()]))
        ttk.Label(self.scope_box, textvariable=self.scope_help, wraplength=325).pack(fill=tk.X, pady=4)
        self.model_box = ttk.Frame(task_box)
        self.model_box.pack(fill=tk.X)
        ttk.Label(self.model_box, text="Model source").pack(anchor="w")
        self.model_selector = ttk.Combobox(self.model_box, textvariable=self.model_source,
                                          values=tuple(MODEL_SOURCES), state="readonly")
        self.model_selector.pack(fill=tk.X)
        self.choose_model_button = ttk.Button(self.model_box, text="Choose local model...", command=self.choose_local_model)
        self.choose_model_button.pack(fill=tk.X, pady=(4, 0))
        self.custom_box = ttk.Frame(self.model_box)
        self._entry(self.custom_box, "model", "Model file (.joblib)")
        ttk.Button(self.custom_box, text="Read model target and scope", command=self.inspect_local_model).pack(fill=tk.X, pady=3)
        ttk.Label(self.custom_box, textvariable=self.local_model_details, wraplength=325).pack(fill=tk.X, pady=3)
        self.mutation_box = ttk.Frame(task_box)
        ttk.Label(self.mutation_box, text="One row per variant; join substitutions with commas.\nExample: A:G1W,B:G1A. Optional ID + Tab + substitutions.", wraplength=325).pack(fill=tk.X, pady=(4, 2))
        self.variants = tk.Text(self.mutation_box, height=4, width=38, wrap=tk.NONE, undo=True)
        self.variants.pack(fill=tk.X)
        import_options = ttk.Frame(self.mutation_box)
        import_options.pack(fill=tk.X, pady=(3, 0))
        self.import_mode = tk.StringVar(value="append")
        ttk.Label(import_options, text="CSV import:").pack(side=tk.LEFT)
        for label, mode in (("Append rows", "append"), ("Replace draft", "replace")):
            ttk.Radiobutton(import_options, text=label, variable=self.import_mode, value=mode).pack(side=tk.LEFT, padx=3)
        ttk.Button(self.mutation_box, text="Import variant CSV...", command=self.import_variants).pack(fill=tk.X, pady=3)
        self._mutation_builder(self.mutation_box)
        self._entry(self.mutation_box, "scope_note", "Why this model fits your system (optional)")
        self.layout = tk.BooleanVar(value=False)
        ttk.Checkbutton(inputs, text="Also generate 2D interface atlases", variable=self.layout).pack(anchor="w", pady=(8, 2))
        ttk.Label(inputs, text="The 3D interface is included with each prediction. Atlases can also be generated from a saved result.", wraplength=340).pack(fill=tk.X, pady=(0, 6))

        output_box = ttk.LabelFrame(inputs, text="3. Save results", padding=6)
        output_box.pack(fill=tk.X, pady=(0, 8))
        self._entry(output_box, "output", "Results folder (leave blank to create one)",
                    str(Path(selected("var_output_dir")) / "affinity") if selected("var_output_dir") else "", browse="directory")
        self.retry = tk.BooleanVar(value=True)
        ttk.Checkbutton(output_box, text="Keep previous runs in separate directories", variable=self.retry).pack(anchor="w", pady=4)
        thermo = ttk.LabelFrame(advanced, text="Thermodynamic conversions", padding=6)
        thermo.pack(fill=tk.X, pady=(0, 8))
        self._entry(thermo, "temperature", "Reference temperature (K)", "298.15")
        self._entry(thermo, "wt_kd", "Measured WT Kd (M), for mutation tasks")
        ttk.Label(thermo, text="Optional WT Kd anchors the mutant Kd estimate; for example, 1e-8 M = 10 nM.", wraplength=325).pack(fill=tk.X, pady=4)
        compute = ttk.LabelFrame(advanced, text="Compute and resources", padding=6)
        compute.pack(fill=tk.X, pady=(0, 8))
        for key, label, value, browse in [
            ("device", "Compute device (auto / cpu / cuda)", "auto", None),
            ("resource_cache", "Resource cache (blank: default)", "", "directory"),
            ("sequence_dir", "Embedding directory override", "", "directory"),
        ]:
            self._entry(compute, key, label, value, browse=browse)
        ttk.Label(compute, textvariable=self.resource_details, wraplength=325).pack(fill=tk.X, pady=6)
        self.repair_button = ttk.Button(compute, text="Manage prediction dependencies...", command=self.repair_dependencies)
        self.repair_button.pack(fill=tk.X, pady=4)
        prior = ttk.LabelFrame(advanced, text="Additional mutation predictor", padding=6)
        prior.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(prior, text="Some mutation models use an additional predictor. Use the matching resources specified by your model provider.", wraplength=325).pack(fill=tk.X)
        for key, label, browse in [("prior_csv", "Predictor CSV (matching variant IDs)", "file"),
                                   ("prior_repository", "Predictor source folder", "directory"),
                                   ("prior_checkpoint", "Predictor checkpoint", "file")]:
            self._entry(prior, key, label, browse=browse)
        other = ttk.LabelFrame(advanced, text="Structure metadata and atlas settings", padding=6)
        other.pack(fill=tk.X)
        for key, label, value, browse in [
            ("source", "Structure source / numbering", "user", None),
            ("complex_id", "Complex ID", "", None),
            ("optcuts_bin", "OptCuts binary", selected("entry_optcuts_bin", "OptCuts_bin"), "file"),
            ("layout_timeout", "Layout timeout per partner (s)", "600", None),
        ]:
            self._entry(other, key, label, value, browse=browse)

        activity = ttk.Frame(self.tabs, padding=6)
        self.tabs.add(activity, text="Activity")
        self.activity_tab = activity
        ttk.Label(activity, text="Progress and troubleshooting details", wraplength=325).pack(fill=tk.X, pady=(0, 5))
        self.log_text = scrolledtext.ScrolledText(activity, width=38, height=15, wrap=tk.WORD, state=tk.DISABLED)
        self.log_text.pack(fill=tk.BOTH, expand=True)
        log_actions = ttk.Frame(activity)
        log_actions.pack(fill=tk.X, pady=5)
        ttk.Button(log_actions, text="Copy diagnostics", command=self.copy_diagnostics).pack(side=tk.LEFT)
        ttk.Button(log_actions, text="Save log...", command=self.save_diagnostics).pack(side=tk.RIGHT)

        resources = ttk.LabelFrame(left, text="4. Start prediction", padding=6)
        resources.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(resources, textvariable=self.model_status, wraplength=340).pack(fill=tk.X, pady=(0, 4))
        self.run_button = ttk.Button(resources, text="Predict", style="Primary.TButton", command=self.start)
        self.run_button.pack(fill=tk.X)
        buttons = ttk.Frame(resources)
        buttons.pack(fill=tk.X, pady=(5, 0))
        self.prepare_button = ttk.Button(buttons, text="Prepare model resources", command=self.prepare_models)
        self.prepare_button.pack(side=tk.LEFT)
        self.cancel_button = ttk.Button(buttons, text="Cancel", command=self.cancel_run, state=tk.DISABLED)
        self.progress = ttk.Progressbar(left, mode="determinate", maximum=100, value=0)
        self.status_label = ttk.Label(left, textvariable=self.status, wraplength=350)
        self.status_label.pack(fill=tk.X, pady=6)

        heading = ttk.Frame(right)
        heading.pack(fill=tk.X, pady=(0, 6))
        self.heading_font = tkfont.nametofont("TkDefaultFont").copy()
        self.heading_font.configure(size=12, weight="bold")
        ttk.Label(heading, text="Prediction results", font=self.heading_font).pack(side=tk.LEFT)
        self.open_button = ttk.Button(heading, text="Open saved result...", command=self.open_result)
        self.open_button.pack(side=tk.RIGHT)
        summary_label = ttk.Label(right, textvariable=self.summary, wraplength=720)
        summary_label.pack(fill=tk.X)
        table_frame = ttk.Frame(right)
        table_frame.pack(fill=tk.X, pady=6)
        self.table = ttk.Treeview(table_frame, columns=("variant", "substitutions", "prediction", "status"), show="headings", height=5)
        for key, title, col_width in [("variant", "Variant", 105), ("substitutions", "Substitutions", 170),
                                      ("prediction", "Prediction", 120), ("status", "Status", 105)]:
            self.table.heading(key, text=title)
            self.table.column(key, width=col_width, minwidth=70, stretch=True)
        self.table.grid(row=0, column=0, sticky="nsew")
        yscroll = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.table.yview)
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=self.table.xview)
        xscroll.grid(row=1, column=0, sticky="ew")
        self.table.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)
        table_frame.columnconfigure(0, weight=1)
        self.table.bind("<<TreeviewSelect>>", self._table_selected)
        selectors = ttk.Frame(right)
        selectors.pack(fill=tk.X)
        ttk.Label(selectors, text="Partner").grid(row=0, column=0, sticky="w", padx=3)
        self.partner_selector = ttk.Combobox(selectors, textvariable=self.partner_label,
            values=self.partner_labels, state="readonly", width=10)
        self.partner_selector.grid(row=1, column=0, sticky="ew", padx=3)
        self.partner_selector.bind("<<ComboboxSelected>>", self._choose_partner)
        self.view_selectors = [self.partner_selector]
        for column, (label, variable, choices, field, values) in enumerate([
            ("View", self.view_label, VIEWS, self.mode, VIEWS),
            ("Colour by", self.descriptor_label, DESCRIPTORS, self.descriptor, DESCRIPTORS)], start=1):
            ttk.Label(selectors, text=label).grid(row=0, column=column, sticky="w", padx=3)
            selector = ttk.Combobox(selectors, textvariable=variable, values=tuple(choices), state="readonly", width=18)
            selector.grid(row=1, column=column, sticky="ew", padx=3)
            selector.bind("<<ComboboxSelected>>", lambda _e, v=variable, f=field, m=values: self._choose_view(v, f, m))
            self.view_selectors.append(selector)
        for column in range(3):
            selectors.columnconfigure(column, weight=1)
        review_actions = ttk.Frame(right)
        review_actions.pack(fill=tk.X, pady=4)
        for index, (text, command) in enumerate([("Export predictions...", self.export_predictions), ("Open result folder", self.open_folder),
                                                ("Reuse inputs", self.reuse_inputs)]):
            button = ttk.Button(review_actions, text=text, command=command)
            button.grid(row=index // 2, column=index % 2, sticky="ew", padx=2, pady=2)
            self.review_buttons.append(button)
        self.layout_button = ttk.Button(review_actions, text="Generate 2D atlases", command=self.retry_layout)
        self.layout_button.grid(row=1, column=1, sticky="ew", padx=2, pady=2)
        self.linked_button = ttk.Button(review_actions, text="Detailed interface review...", command=self.open_linked_review)
        self.linked_button.grid(row=2, column=0, columnspan=2, sticky="ew", padx=2, pady=2)
        review_actions.columnconfigure((0, 1), weight=1)
        self.show_main_button = ttk.Button(review_actions, text="Show in main window", command=self.show_in_main)
        if host is not None and hasattr(host, "receive_affinity_result"):
            self.show_main_button.grid(row=3, column=0, columnspan=2, sticky="ew", padx=2, pady=2)
        error_label = ttk.Label(right, textvariable=self.errors, foreground="#a12222", wraplength=720)
        details_label = ttk.Label(right, textvariable=self.details, wraplength=720)
        self.show_details = tk.BooleanVar(value=False)
        def toggle_details():
            if self.show_details.get():
                details_label.pack(fill=tk.X, pady=4, before=self.plot)
            else:
                details_label.pack_forget()
        details_toggle = ttk.Checkbutton(right, text="Show result details", variable=self.show_details,
                                         command=toggle_details)
        details_toggle.pack(anchor=tk.W, pady=2)
        def show_errors(*_args):
            if self.errors.get().strip():
                error_label.pack(fill=tk.X, pady=3, before=details_toggle)
            else:
                error_label.pack_forget()
        self.errors.trace_add("write", show_errors)
        def fit_results(event):
            # Keep result controls and the canvas usable in compact windows.
            self.table.configure(height=1 if event.height < 560 else 3 if event.height < 720 else 5)
            for label in (summary_label, error_label, details_label):
                label.configure(wraplength=max(200, event.width - 16))
            buttons = [*self.review_buttons, self.layout_button, self.linked_button]
            if host is not None and hasattr(host, "receive_affinity_result"):
                buttons.append(self.show_main_button)
            columns = 3 if event.width >= 640 else 2
            for column in range(3):
                review_actions.columnconfigure(column, weight=1 if column < columns else 0)
            for index, button in enumerate(buttons):
                button.grid(row=index // columns, column=index % columns, columnspan=1, sticky="ew", padx=2, pady=2)
        right.bind("<Configure>", fit_results)
        self.plot = ttk.Frame(right)
        self.plot.pack(fill=tk.BOTH, expand=True)
        self._poll_id = self.window.after(100, self._poll)
        for variable in (self.task, self.mutation_context, self.model_source):
            variable.trace_add("write", lambda *_: self._task_changed())
        for key in ("model", "resource_cache", "prior_csv", "prior_repository"):
            self.vars[key].trace_add("write", lambda *_: self._task_changed())
        self.vars["structure"].trace_add("write", lambda *_: self._invalidate_structure())
        for key in ("structure", "chains_a", "chains_b"):
            self.vars[key].trace_add("write", lambda *_: self._refresh_main_input_context())
        self.window.bind("<FocusIn>", lambda _event: self._refresh_main_input_context(), add="+")
        self._refresh_main_input_context()
        self._task_changed()
        self._review_state()
        if analysis:
            self.load(analysis)
        elif self.vars["structure"].get().strip():
            self.read_structure()

    def _main_window_values(self):
        fields = {"structure": "var_input_path", "chains_a": "var_chain_a", "chains_b": "var_chain_b",
                  "output": "var_output_dir"}
        return {key: str(getattr(self.host, name).get()).strip() if self.host is not None and hasattr(self.host, name) else ""
                for key, name in fields.items()}

    def _refresh_main_input_context(self):
        values = self._main_window_values()
        available = bool(values["structure"])
        self.use_main_input_button.configure(state=tk.NORMAL if available and not self._working and not self.closing else tk.DISABLED)
        if not available:
            self.main_input_context.set("Select a structure in the main window to reuse its input here.")
            return
        description = f"{Path(values['structure']).name} | {values['chains_a'] or '?'} / {values['chains_b'] or '?'}"
        matches = all(values[key] == self.vars[key].get().strip() for key in ("structure", "chains_a", "chains_b"))
        self.main_input_context.set("Main window: " + description +
            ("\nStructure and partners match this setup." if matches else "\nUse its structure and partners for your next prediction."))

    def use_main_window_input(self):
        if self._working or self.closing:
            return
        values = self._main_window_values()
        if not values["structure"]:
            self._refresh_main_input_context()
            self.status.set("Choose a structure in the main window, then use its input here.")
            return
        from topoppi_affinity.workflow import inspect_structure

        try:
            records = inspect_structure(values["structure"])
        except (OSError, ValueError, KeyError) as exc:
            self.errors.set(f"Check the structure selected in the main window: {exc}")
            return
        for key in ("structure", "chains_a", "chains_b"):
            self.vars[key].set(values[key])
        if values["output"]:
            self.vars["output"].set(str(Path(values["output"]) / "affinity"))
        self.errors.set("")
        self.read_structure(records=records)
        self._refresh_main_input_context()
        self.tabs.select(0)
        message = f"Using {Path(values['structure']).name} from the main window."
        if self.variants.get("1.0", "end-1c").strip():
            message += " Variant rows are retained; review their residue identities for this structure."
        if self.document is not None:
            message += " The displayed result is unchanged."
        self.status.set(message)

    def open_linked_review(self):
        if self.document is None or self.document.interface is None:
            self.errors.set("Open an analysis with a saved native interface first.")
            return
        from topoppi_affinity.atlas_review.gui import ReviewWindow
        from topoppi_affinity.atlas_review.session import ReviewSession
        has_atlas = any(item.get("status") == "ok" for item in self.document.report.get("layout", {}).values())
        session = ReviewSession(self.document, condition="A" if has_atlas else "G")
        selection = self.table.selection()
        if selection:
            session.select_variant(selection[0])
        session.side = int(self.side.get())
        ReviewWindow(self.window, session)

    def _scroll_tab(self, tabs, title):
        from topoppi.gui_app.scrolling import reveal_focused_widget, scroll_settings_page

        outer = ttk.Frame(tabs)
        tabs.add(outer, text=title)
        canvas = tk.Canvas(outer, highlightthickness=0, width=355)
        scrollbar = ttk.Scrollbar(outer, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        content = ttk.Frame(canvas, padding=6)
        item = canvas.create_window(0, 0, anchor="nw", window=content)
        content.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(item, width=e.width))
        # Bind only this panel's descendants, so scrolling never captures another window.
        def wheel(event):
            if event.num in (4, 5):
                amount = -1 if event.num == 4 else 1
            else:
                amount = -1 if event.delta > 0 else 1
            widget = event.widget
            if str(widget).startswith(str(content)) and not isinstance(widget, (tk.Text, tk.Listbox, ttk.Combobox)):
                canvas.yview_scroll(amount * 3, "units")
        self.window.bind("<MouseWheel>", wheel, add="+")
        self.window.bind("<Button-4>", wheel, add="+")
        self.window.bind("<Button-5>", wheel, add="+")
        self.window.bind("<FocusIn>", lambda event: reveal_focused_widget(canvas, content, event.widget), add="+")
        for sequence, direction in (("<Prior>", -1), ("<Next>", 1)):
            self.window.bind(sequence, lambda event, direction=direction: scroll_settings_page(
                canvas, content, event, direction), add="+")
        return content

    def _entry(self, frame, key, label, value="", browse=None):
        ttk.Label(frame, text=label).pack(anchor="w", pady=(4, 0))
        row = ttk.Frame(frame)
        row.pack(fill=tk.X)
        variable = self.vars[key] = tk.StringVar(value=value)
        ttk.Entry(row, textvariable=variable, width=40).pack(side=tk.LEFT, fill=tk.X, expand=True)
        if browse:
            def choose():
                path = (filedialog.askdirectory(parent=self.window) if browse == "directory"
                        else filedialog.askopenfilename(parent=self.window))
                if path:
                    variable.set(path)
                    if key == "structure":
                        self.read_structure()
            ttk.Button(row, text="...", width=3, command=choose).pack(side=tk.LEFT)

    def _resource_plan(self, task, *, cache_dir=None, prior_override=False):
        from topoppi_affinity.resources import resource_status

        status = resource_status(cache_dir=cache_dir)
        model = status["models"].get(task, {})
        missing = list(model.get("missing_resources", []))
        if task is None:
            missing = []
        elif prior_override:
            missing = ([] if model["installed"] else [model["name"] + " model"])
            if not status["esmc"]["installed"]:
                missing.append("ESMC-600M sequence weights")
        return status, missing

    def _workflow_task(self):
        """Translate form choices to the unchanged scientific workflow contract."""
        if MODEL_SOURCES[self.model_source.get()] == "local":
            return None
        return "absolute" if TASKS[self.task.get()] == "absolute" else (self.mutation_context.get() or "mutation")

    def _needs_local_model(self):
        task = self._workflow_task()
        if task is None:
            return not self.vars["model"].get().strip()
        model = (self.resource_state or {}).get("models", {}).get(task, {})
        return bool(model and not model.get("installed") and not model.get("url"))

    def _read_local_model(self):

        from topoppi_affinity.inference import model_info

        path = self.vars["model"].get().strip()
        if not path or not Path(path).expanduser().is_file():
            raise ValueError("Choose an existing model.json or .joblib file in step 2")
        try:
            from topoppi_affinity.production import load_model
            info = model_info(load_model(Path(path).expanduser()))
        except Exception as exc:
            raise ValueError(f"Cannot read this local model: {exc}. Choose a compatible TopoPPI model.") from exc
        if info.get("target") not in {"pkd", "ddg"}:
            raise ValueError("The local model must declare a pKd or ddG prediction target")
        self._local_model_path, self._local_model_info = path, info
        return info

    def inspect_local_model(self):
        """An explicit selection adopts the artifact's target, retaining variant edits."""
        if self._working or self.closing:
            return
        try:
            info = self._read_local_model()
        except ValueError as exc:
            self._local_model_info = None
            self.errors.set(str(exc))
            self._task_changed()
            return
        self.task.set("Binding affinity (pKd)" if info["target"] == "pkd" else "Mutation effects (ddG)")
        if info["kind"] in {"known", "unseen"}:
            self.mutation_context.set(info["kind"])
        self.model_source.set("Local model file")
        self.errors.set("")
        self._task_changed()
        retained = " Variant rows are retained; pKd prediction uses the unmodified complex." if info["target"] == "pkd" and self.variants.get("1.0", "end-1c").strip() else ""
        self.status.set(f"Model target: {'pKd' if info['target'] == 'pkd' else 'ddG'}. Review the structure and inputs, then select Predict." + retained)

    def choose_local_model(self):
        """Select a supplied fitted model while retaining the current inputs."""
        if self._working or self.closing:
            return
        path = filedialog.askopenfilename(parent=self.window, title="Choose a binding prediction model",
                                          filetypes=[("Fitted model", "*.joblib *.json"), ("All files", "*")])
        if path:
            self.vars["model"].set(path)
            self.model_source.set("Local model file")
            self.inspect_local_model()
            self.tabs.select(0)

    def _task_changed(self):
        prediction_task = TASKS[self.task.get()]
        task = self._workflow_task()
        self.task_help.set(TASK_HELP[prediction_task])
        self.scope_help.set(SCOPE_HELP[self.mutation_context.get()])
        self.mutation_model_choice.set(next(label for label, value in MUTATION_MODELS.items()
                                           if value == self.mutation_context.get()))
        self.scope_selector.configure(state="readonly" if prediction_task == "mutation" else tk.DISABLED)
        if self.vars["model"].get().strip() != self._local_model_path:
            self._local_model_info = None
        info = self._local_model_info
        self.custom_box.pack_forget()
        self.mutation_box.pack_forget()
        if task is None:
            self.custom_box.pack(fill=tk.X)
        if prediction_task == "mutation":
            self.mutation_box.pack(fill=tk.X)
        if task is None:
            if info:
                target = "pKd" if info["target"] == "pkd" else "ddG (kcal/mol)"
                scope = info["scope"]
                if info["target"] == "ddg" and info["kind"] not in {"known", "unseen"}:
                    scope += ". Known/unseen applicability is not declared by this artifact; your selection records your assessment."
                self.local_model_details.set(f"Target: {target}. {scope}")
            else:
                self.local_model_details.set("Read the model target and scope before predicting. Your variant rows will be retained.")
        if task == "":
            self.resource_state = None
            self.model_status.set("Choose known or unseen model applicability before preparing a mutation prediction.")
            self.resource_details.set("Model resources depend on the selected applicability.")
            self.run_button.configure(text="Predict")
            self.prepare_button.pack_forget()
            return
        cache_dir = self.vars["resource_cache"].get().strip() or None
        try:
            status, missing = self._resource_plan(task, cache_dir=cache_dir,
                prior_override=bool(self.vars["prior_csv"].get().strip() or self.vars["prior_repository"].get().strip()))
            self.resource_state = status
            detail = f"Cache: {status['cache_dir']}"
            if task is not None:
                model = status["models"][task]
                detail += f"\n{model['scope']}"
                if model.get("prior_license"):
                    detail += f"\nPrior terms: {model['prior_license']}"
            else:
                detail += "\nLocal models use their supplied training scope and prior configuration."
            self.resource_details.set(detail)
            if missing:
                size = (f" ESMC download: {status['esmc']['size_bytes'] / 1e9:.2f} GB."
                        if not status['esmc']['installed'] else "")
                self.model_status.set("To prepare: " + ", ".join(missing) + "." + size)
            else:
                self.model_status.set("Resources ready. The default physical models run locally on CPU."
                                      if task in {"absolute", "mutation"} else "Prepared model resources are ready.")
            self.run_button.configure(text="Prepare and predict" if missing else "Predict")
            if task is None:
                self.model_status.set("Local model resources are checked when prediction starts." if self.vars["model"].get().strip()
                                      else "Choose a local model file to enable prediction.")
                self.run_button.configure(text="Choose local model..." if self._needs_local_model() else "Predict")
                self.prepare_button.pack_forget()
                expected = "pkd" if prediction_task == "absolute" else "ddg"
                if info and info["target"] != expected:
                    self.model_status.set(f"Model target does not match this task. Choose a {expected} model, or change the prediction task.")
                elif info and prediction_task == "mutation" and info["kind"] in {"known", "unseen"} and info["kind"] != self.mutation_context.get():
                    self.model_status.set(f"This model is for {info['kind']} systems. Match the applicability selection or choose another model.")
            elif self._needs_local_model():
                self.model_status.set("Choose a local .joblib model to enable prediction. Saved results are ready to open.")
                self.run_button.configure(text="Choose local model...")
                self.prepare_button.pack_forget()
            else:
                if task in {"absolute", "mutation"}:
                    self.prepare_button.pack_forget()
                else:
                    self.prepare_button.pack(fill=tk.X, pady=(4, 0))
                dependencies = self._missing_dependencies(task)
                if dependencies:
                    self.model_status.set("Install prediction dependencies: " + ", ".join(dependencies) + ".")
                    self.run_button.configure(text="Set up prediction...")
        except (OSError, ValueError, KeyError) as exc:
            self.resource_state = None
            self.model_status.set("Check the resource cache in Advanced.")
            self.resource_details.set(str(exc))
            self.run_button.configure(text="Predict")

    def _review_state(self):
        available = self.document is not None
        has_interface = available and getattr(self.document, "interface", None) is not None
        active = self._working
        for button in self.review_buttons:
            button.configure(state=tk.NORMAL if available else tk.DISABLED)
        self.review_buttons[2].configure(state=tk.NORMAL if available and not active else tk.DISABLED)
        self.layout_button.configure(state=tk.NORMAL if has_interface and not active else tk.DISABLED)
        self.linked_button.configure(state=tk.NORMAL if has_interface else tk.DISABLED)
        self.show_main_button.configure(state=tk.NORMAL if available else tk.DISABLED)
        for selector in self.view_selectors:
            selector.configure(state="readonly" if has_interface else tk.DISABLED)

    def _choose_view(self, label, field, mapping):
        field.set(mapping[label.get()])
        self.draw()

    def _mutation_builder(self, parent):
        builder = ttk.LabelFrame(parent, text="Add a specific substitution", padding=4)
        builder.pack(fill=tk.X, pady=4)
        self.mutation_vars = {key: tk.StringVar(value="") for key in ("chain", "position", "insertion", "wt", "mt")}
        self.residue_choice = tk.StringVar()
        top = ttk.Frame(builder)
        top.pack(fill=tk.X)
        ttk.Label(top, text="Chain").pack(side=tk.LEFT)
        self.mutation_chain = ttk.Combobox(top, textvariable=self.mutation_vars["chain"], width=7, state="readonly")
        self.mutation_chain.pack(side=tk.LEFT, padx=3)
        self.mutation_chain.bind("<<ComboboxSelected>>", lambda _e: self._residues_for_chain())
        self.mutation_residue = ttk.Combobox(top, textvariable=self.residue_choice, width=20, state="readonly")
        self.mutation_residue.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.mutation_residue.bind("<<ComboboxSelected>>", lambda _e: self._select_residue())
        fields = ttk.Frame(builder)
        fields.pack(fill=tk.X)
        for key, label, width in [("position", "Residue", 7), ("insertion", "Ins. code", 7), ("wt", "From", 4), ("mt", "To", 4)]:
            box = ttk.Frame(fields)
            box.pack(side=tk.LEFT, padx=2)
            ttk.Label(box, text=label).pack()
            if key == "mt":
                ttk.Combobox(box, textvariable=self.mutation_vars[key], values=AMINO_ACIDS, state="readonly", width=width).pack()
            else:
                ttk.Entry(box, textvariable=self.mutation_vars[key], width=width,
                          state="readonly" if key == "wt" else "normal").pack()
        self.mutation_vars["position"].trace_add("write", lambda *_: self._lookup_wt())
        self.mutation_vars["insertion"].trace_add("write", lambda *_: self._lookup_wt())
        row = ttk.Frame(builder)
        row.pack(fill=tk.X, pady=3)
        ttk.Button(row, text="New variant row", command=self.add_mutation).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(row, text="Add to selected row", command=lambda: self.add_mutation(joint=True)).pack(side=tk.LEFT, expand=True, fill=tk.X)

    def read_structure(self, records=None):
        from topoppi_affinity.workflow import inspect_structure
        try:
            self.structure_records = records if records is not None else inspect_structure(self.vars["structure"].get().strip())
            self.chain_choices.delete(0, tk.END)
            for chain in self.structure_records["chains"]:
                self.chain_choices.insert(tk.END, f"{chain['id']}  ({chain['residue_count']} residues)")
            chains = [chain["id"] for chain in self.structure_records["chains"]]
            self.chain_status.set(f"{len(chains)} protein chains found. Select chains below and assign each binding partner.")
            selected_groups = [set(self.vars[f"chains_{side}"].get().replace(" ", "").split(",")) for side in ("a", "b")]
            if len(chains) == 2 and (not all(group and group <= set(chains) for group in selected_groups)
                                    or selected_groups[0] & selected_groups[1]):
                self.vars["chains_a"].set(chains[0])
                self.vars["chains_b"].set(chains[1])
            self.mutation_chain.configure(values=chains)
            if chains:
                partner_chains = (chain.strip() for side in ("a", "b")
                                  for chain in self.vars[f"chains_{side}"].get().split(","))
                self.mutation_vars["chain"].set(next((chain for chain in partner_chains if chain in chains), chains[0]))
                self._residues_for_chain()
            if self.errors.get().startswith("Cannot read structure:"):
                self.errors.set("")
            self.status.set("Structure loaded. Check the two partner groups, then choose a prediction task.")
        except (OSError, ValueError, KeyError) as exc:
            self.structure_records = None
            self.chain_choices.delete(0, tk.END)
            self.chain_status.set("Choose a readable PDB or mmCIF file, then select Read chains.")
            self.errors.set(f"Cannot read structure: {exc}")

    def _invalidate_structure(self):
        self.chain_status.set("Select Read chains to inspect this structure and assign its binding partners.")
        self.structure_records = None
        self.chain_choices.delete(0, tk.END)
        self.mutation_chain.configure(values=())
        self.mutation_residue.configure(values=())
        self.residue_choice.set("")
        for key in ("chain", "position", "insertion", "wt"):
            self.mutation_vars[key].set("")

    def _assign_chains(self, side):
        if self.structure_records is None:
            self.read_structure()
        if self.structure_records is not None:
            ids = [self.structure_records["chains"][i]["id"] for i in self.chain_choices.curselection()]
            if ids:
                self.vars[f"chains_{side}"].set(",".join(ids))
                other = "b" if side == "a" else "a"
                remaining = [chain.strip() for chain in self.vars[f"chains_{other}"].get().split(",")
                             if chain.strip() and chain.strip() not in ids]
                self.vars[f"chains_{other}"].set(",".join(remaining))
                self.status.set(f"Partner {side.upper()}: {', '.join(ids)}. Check both partner groups before prediction.")

    def _residues_for_chain(self):
        chain = self.mutation_vars["chain"].get()
        residues = (self.structure_records or {}).get("residues", [])
        self._residue_options = [r for r in residues if r["chain"] == chain]
        labels = [f"{r['position']}{r['insertion']} | {r['wt']}" for r in self._residue_options]
        self.mutation_residue.configure(values=labels)
        self.residue_choice.set(labels[0] if labels else "")
        self._select_residue()

    def _select_residue(self):
        index = self.mutation_residue.current()
        if index >= 0:
            residue = self._residue_options[index]
            for key in ("position", "insertion", "wt"):
                self.mutation_vars[key].set(str(residue[key]))
        else:
            for key in ("position", "insertion", "wt"):
                self.mutation_vars[key].set("")

    def _lookup_wt(self):
        values = {key: variable.get().strip() for key, variable in self.mutation_vars.items()}
        residues = (self.structure_records or {}).get("residues", [])
        found = next((r for r in residues if r["chain"] == values["chain"] and
                      str(r["position"]) == values["position"] and r["insertion"] == values["insertion"]), None)
        self.mutation_vars["wt"].set(found["wt"] if found else "")

    def add_mutation(self, joint=False):
        values = {key: variable.get().strip() for key, variable in self.mutation_vars.items()}
        if not all(values[key] for key in ("chain", "position", "wt", "mt")):
            self.errors.set("Read a structure, select an existing residue and choose its replacement amino acid.")
            return
        if values["wt"] == values["mt"]:
            self.errors.set("Choose a replacement amino acid different from the structure residue.")
            return
        mutation = f"{values['chain']}:{values['wt']}{values['position']}{values['insertion']}{values['mt']}"
        if joint and self.variants.get("1.0", "end-1c").strip():
            line = self.variants.index(tk.INSERT).split(".")[0]
            content = self.variants.get(f"{line}.0", f"{line}.end")
            if content.strip():
                self.variants.insert(f"{line}.end", "," + mutation)
            else:
                self.variants.insert(f"{line}.0", mutation)
        else:
            content = self.variants.get("1.0", "end-1c")
            self.variants.insert(tk.END, ("\n" if content and not content.endswith("\n") else "") + mutation)
        self.errors.set("")

    def import_variants(self):
        path = filedialog.askopenfilename(parent=self.window, title="Variant CSV: variant_id, mutations",
                                          filetypes=[("CSV", "*.csv"), ("All files", "*")])
        if not path:
            return
        from topoppi_affinity.workflow import variants_table
        try:
            frame = variants_table(variants=path)
            existing_text = self.variants.get("1.0", "end-1c")
            append = self.import_mode.get() == "append"
            if append and existing_text.strip():
                existing = variant_lines(existing_text)
                duplicates = set(existing.variant_id) & set(frame.variant_id)
                if duplicates:
                    raise ValueError("These IDs already exist in the draft: " + ", ".join(sorted(duplicates)) +
                                     ". Edit the conflicting IDs, or choose Replace draft to import this file on its own.")
            incoming = "\n".join(f"{r.variant_id}\t{r.mutations}" for r in frame.itertuples())
            self.variants.edit_separator()
            self.variants.configure(autoseparators=False)
            if append:
                separator = "\n" if existing_text and not existing_text.endswith("\n") else ""
                self.variants.insert(tk.END, separator + incoming)
            else:
                self.variants.delete("1.0", tk.END)
                self.variants.insert("1.0", incoming)
            self.variants.configure(autoseparators=True)
            self.variants.edit_separator()
            if self.errors.get().startswith("Cannot import variants:"):
                self.errors.set("")
            action = "Appended" if append else "Replaced draft with"
            self.status.set(f"{action} {len(frame)} variants; sample IDs are preserved. Ctrl+Z undoes this import.")
        except (OSError, ValueError) as exc:
            self.errors.set(f"Cannot import variants: {exc}")

    def _values(self):
        values = {key: variable.get().strip() for key, variable in self.vars.items()}
        if not all(values[key] for key in ("structure", "chains_a", "chains_b")):
            raise ValueError("Structure and both complete partner groups are required")
        prediction_task = TASKS[self.task.get()]
        values["task"] = self._workflow_task()
        if values["task"] is None:
            info = self._read_local_model()
            expected = "pkd" if prediction_task == "absolute" else "ddg"
            if info["target"] != expected:
                raise ValueError(f"The local model predicts {info['target']}, but this task needs {expected}. Change the task or choose a matching model.")
            if prediction_task == "mutation" and info["kind"] in {"known", "unseen"} and info["kind"] != self.mutation_context.get():
                raise ValueError(f"The local model is for {info['kind']} systems. Change Mutation model applicability or choose a matching model.")
        else:
            values["model"] = ""
        if prediction_task == "absolute":
            for key in ("wt_kd", "prior_csv", "prior_repository", "prior_checkpoint", "scope_note"):
                values[key] = ""
        if not self.layout.get():
            values["layout_timeout"] = "600"
        for key, label, unit, example in (("temperature", "Reference temperature", "kelvin", "310"),
                                           ("layout_timeout", "Layout timeout", "seconds", "600")):
            try:
                values[key] = float(values[key])
            except ValueError as exc:
                raise ValueError(f"{label} must be a number in {unit}, for example {example}. Enter the number without a unit suffix.") from exc
            if not math.isfinite(values[key]) or values[key] <= 0:
                raise ValueError(f"{label} must be finite and positive")
        try:
            values["wt_kd"] = float(values["wt_kd"]) if values["wt_kd"] else None
        except ValueError as exc:
            raise ValueError("Measured WT Kd must be a number in M, for example 1e-8. Enter the number without a unit suffix.") from exc
        if values["wt_kd"] is not None and (not math.isfinite(values["wt_kd"]) or values["wt_kd"] <= 0):
            raise ValueError("Measured WT Kd must be finite and positive (M)")
        for key in ("model", "output", "sequence_dir", "prior_csv", "prior_repository", "prior_checkpoint", "complex_id"):
            values[key] = values[key] or None
        if values["output"]:
            output = Path(values["output"]).expanduser()
            if output.exists() and not output.is_dir():
                raise ValueError("Results folder points to a file. Choose a folder for the new analysis.")
            if (output / "analysis.json").exists() and not self.retry.get():
                raise ValueError("This folder already contains an analysis. Choose another results folder or enable "
                                 "'Keep previous runs in separate directories'.")
        self._check_prior_paths(values)
        values["variants"] = (None if prediction_task == "absolute" else variant_lines(self.variants.get("1.0", tk.END)))
        if prediction_task == "mutation" and values["variants"] is None:
            raise ValueError("Add at least one substitution in step 2, or import a variant CSV")
        self._check_structure_inputs(values)
        values["layout"], values["retry"] = self.layout.get(), self.retry.get()
        if prediction_task == "mutation":
            values["mutation_context"] = self.mutation_context.get() or None
        return values

    @staticmethod
    def _check_prior_paths(values):
        if values["prior_csv"] and values["prior_repository"]:
            raise ValueError("Choose a prior CSV or a source repository")
        for key, label in (("prior_csv", "Prior CSV"), ("prior_checkpoint", "Prior checkpoint")):
            if values[key] and not Path(values[key]).expanduser().is_file():
                raise ValueError(f"{label} file does not exist: {values[key]}. Update it in Advanced.")
        if values["prior_repository"] and not Path(values["prior_repository"]).expanduser().is_dir():
            raise ValueError(f"Prior repository directory does not exist: {values['prior_repository']}. Update it in Advanced.")

    def _check_structure_inputs(self, values):
        from topoppi.interface import chain_group

        from topoppi_affinity.features import parse_mutations
        from topoppi_affinity.workflow import inspect_structure

        records = self.structure_records or inspect_structure(values["structure"])
        groups = [chain_group(values[key]) for key in ("chains_a", "chains_b")]
        if set(groups[0]) & set(groups[1]):
            raise ValueError("Assign each chain to one partner; the two partner groups overlap")
        available = {chain["id"] for chain in records["chains"]}
        missing = (set(groups[0]) | set(groups[1])) - available
        if missing:
            raise ValueError("Choose chains present in this structure. Missing: " + ", ".join(sorted(missing)))
        variants = values["variants"]
        if variants is None:
            return
        selected = set(groups[0]) | set(groups[1])
        residues = {(r["chain"], r["position"], r["insertion"]): r["wt"] for r in records["residues"] if r["chain"] in selected}
        errors = []
        for row in variants.itertuples():
            try:
                mutations = parse_mutations(row.mutations)
                if not mutations and (values["task"] in {"known", "unseen"} or TASKS[self.task.get()] == "mutation"):
                    raise ValueError("Add a substitution")
                seen = set()
                for mutation in mutations:
                    key = (mutation["chain"], mutation["position"], mutation["insertion"])
                    if key in seen:
                        raise ValueError("Put alternative substitutions at the same site in separate rows")
                    seen.add(key)
                    if key not in residues:
                        raise ValueError(f"Residue {key} is absent from the selected partner chains")
                    if residues[key] != mutation["wt"]:
                        raise ValueError(f"Residue {key} is {residues[key]} in the structure; update the starting amino acid")
            except (ValueError, KeyError) as exc:
                errors.append(f"{row.variant_id}: {exc}")
        if len(errors) == len(variants):
            raise ValueError("Check the variant rows before prediction. " + "; ".join(errors))

    def _busy(self, active):
        self._working = active
        def controls(widget):
            for child in widget.winfo_children():
                if child is self.activity_tab:
                    continue
                if isinstance(child, (ttk.Entry, ttk.Combobox, ttk.Button, ttk.Checkbutton, ttk.Radiobutton, tk.Listbox, tk.Text)):
                    yield child
                yield from controls(child)
        if active:
            self._input_states = {widget: str(widget.cget("state")) for widget in controls(self.tabs)}
            for widget in self._input_states:
                widget.configure(state=tk.DISABLED)
        else:
            for widget, state in self._input_states.items():
                widget.configure(state=state)
            self._input_states.clear()
        for button in (self.run_button, self.prepare_button, self.layout_button, self.open_button, self.repair_button, self.example_button, self.open_setup_button):
            button.configure(state=tk.DISABLED if active else tk.NORMAL)
        self.cancel_button.configure(state=tk.NORMAL if active else tk.DISABLED)
        if active:
            self.cancel_button.pack(side=tk.RIGHT, padx=(4, 0))
            self.progress.pack(fill=tk.X, pady=(6, 0), before=self.status_label)
            self.progress.configure(mode="indeterminate")
            self.progress.start()
        else:
            self.cancel_button.pack_forget()
            self.progress.stop()
            self.progress.pack_forget()
            self.progress.configure(mode="determinate", maximum=100, value=0)
            self._task_changed()
        self._review_state()
        self._refresh_main_input_context()

    def _launch(self, action):
        if self.worker is not None and self.worker.is_alive():
            return
        self.cancel.clear()
        self._busy(True)
        self._append_activity("Run started")
        def work():
            try:
                action()
            except InterruptedError:
                self.events.put(("cancelled", "Cancelled; completed files have been kept."))
            except Exception as exc:
                self.events.put(("log", traceback.format_exc()))
                self.events.put(("error", f"{type(exc).__name__}: {exc}"))
            finally:
                self.events.put(("finished", None))
        # Release discarded Tk canvas cycles on the UI thread before a worker
        # allocation can trigger their image and variable finalizers.
        gc.collect()
        self.worker = threading.Thread(target=work, name="topoppi-affinity", daemon=False)
        self.worker.start()

    def _can_compute(self):
        installer = getattr(self, "_affinity_installer", None)
        if installer is not None and installer.worker is not None and installer.worker.is_alive():
            self.errors.set("Dependency installation is running. Wait for it to finish, then restart TopoPPI before prediction.")
            return False
        owner = self.host or self
        if getattr(installer, "restart_required", False) or getattr(owner, "_affinity_restart_required", False):
            self.errors.set("Dependencies changed in this session. Use Restart TopoPPI in the installer before prediction.")
            if installer is not None and installer.window.winfo_exists():
                installer.window.lift()
            return False
        return True

    @staticmethod
    def _missing_dependencies(task):
        if task in {None, "absolute", "mutation"}:
            return []
        required = {"torch": "PyTorch", "esm": "ESMC", "transformers": "Transformers"}
        if task == "known":
            required["torch_geometric"] = "PyTorch Geometric"
        return [name for module, name in required.items() if importlib.util.find_spec(module) is None]

    def _prediction_dependencies_ready(self):
        task = self._workflow_task()
        missing = self._missing_dependencies(task)
        if missing:
            from topoppi_affinity.installer import open_installer
            self.status.set("Install prediction dependencies, then restart to continue with these inputs.")
            self.errors.set("Required for this task: " + ", ".join(missing) + ". The setup window will guide installation and restart.")
            open_installer(self, download=True, tasks=(task,))
            return False
        return True

    def start(self):
        if self.worker is not None and self.worker.is_alive():
            return
        if not self._can_compute():
            return
        self._task_changed()
        if self._needs_local_model():
            self.choose_local_model()
            return
        if not self._prediction_dependencies_ready():
            return
        try:
            values = self._values()
        except (ValueError, OSError, KeyError) as exc:
            self.status.set("Review the input error, then predict again.")
            self.errors.set(str(exc))
            messagebox.showerror("Check prediction inputs", str(exc), parent=self.window)
            return
        self.document = None
        self.table.delete(*self.table.get_children())
        self._clear_plot()
        self.summary.set("Preparing a new prediction. Results will appear here as they become available.")
        self.details.set("")
        self.errors.set("")
        self.status.set("Checking resources before prediction...")
        self._review_state()
        values["cache_dir"] = values.pop("resource_cache") or None
        def action():
            from topoppi_affinity.resources import download_models
            from topoppi_affinity.workflow import analyze

            prior_override = bool(values["prior_csv"] or values["prior_repository"])
            _, missing = self._resource_plan(values["task"], cache_dir=values["cache_dir"], prior_override=prior_override)
            if missing:
                options = {"include_priors": False} if prior_override else {}
                download_models(tasks=(values["task"],) if values["task"] else (), include_esmc=True,
                    cache_dir=values["cache_dir"], cancel_event=self.cancel,
                    progress=lambda event: self.events.put(("download", event)), **options)
                if self.cancel.is_set():
                    raise InterruptedError("Resource preparation cancelled")
                self.events.put(("resources_ready", "Resources ready. Starting prediction..."))
            report = analyze(**values, cancel_event=self.cancel, log=lambda text: self.events.put(("log", text)),
                             event_callback=lambda event: self.events.put(("stage", event)))
            self.events.put(("result", Path(report["output_directory"]) / "analysis.json"))
        self._launch(action)

    def repair_dependencies(self):
        from topoppi_affinity.installer import open_installer

        open_installer(self, download=False)

    def cancel_run(self):
        self.cancel.set()
        self.status.set("Cancellation requested. The active computation is finishing; completed files will be kept.")
        self.cancel_button.configure(state=tk.DISABLED)

    def prepare_models(self):
        if not self._can_compute():
            return
        self._task_changed()
        if self._needs_local_model():
            self.choose_local_model()
            return
        if not self._prediction_dependencies_ready():
            return
        task = self._workflow_task()
        if task in {None, "absolute", "mutation"}:
            return
        try:
            self._check_prior_paths({key: self.vars[key].get().strip() if task != "absolute" else ""
                                     for key in ("prior_csv", "prior_repository", "prior_checkpoint")})
        except ValueError as exc:
            self.errors.set(str(exc))
            return
        self.errors.set("")
        self.status.set("Preparing model resources; ESMC is downloaded from its official source.")
        cache_dir = self.vars["resource_cache"].get().strip() or None
        prior_override = bool(self.vars["prior_csv"].get().strip() or self.vars["prior_repository"].get().strip())
        def action():
            from topoppi_affinity.resources import download_models
            options = {"include_priors": False} if prior_override else {}
            download_models(tasks=(task,) if task else (), include_esmc=True, cache_dir=cache_dir, cancel_event=self.cancel,
                            progress=lambda event: self.events.put(("download", event)), **options)
            self.events.put(("resources_ready", "Model resources are ready. Select Predict to run this structure."))
        self._launch(action)

    def load_example(self):
        if not self._can_compute():
            return
        directory = filedialog.askdirectory(parent=self.window,
            title="Choose a folder for the installed 1BVK example", mustexist=False)
        if not directory:
            return
        from topoppi_affinity.resources import export_example

        try:
            example = export_example(directory, reuse=True)
        except (OSError, ValueError) as exc:
            self.errors.set(f"Cannot load the installed example: {exc}")
            self.status.set("Example not loaded. Your current inputs and existing files are retained.")
            return
        for key, value in {"structure": example["structure"], "chains_a": example["chains_a"],
                           "chains_b": example["chains_b"], "output": str(Path(directory) / "results"),
                           "model": "", "wt_kd": "", "temperature": "298.15", "device": "auto",
                           "prior_csv": "", "prior_repository": "", "prior_checkpoint": "", "sequence_dir": "",
                           "source": "software_example", "complex_id": "1BVK",
                           "scope_note": "1BVK antibody-lysozyme software example"}.items():
            self.vars[key].set(value)
        self.task.set("Binding affinity (pKd)")
        self.model_source.set("Prepared task model")
        self.variants.delete("1.0", tk.END)
        self.layout.set(False)
        self.retry.set(True)
        self.mode.set("surface")
        self.side.set("0")
        self.descriptor.set("none")
        self.document = None
        self.table.delete(*self.table.get_children())
        self._clear_plot()
        self.summary.set("1BVK antibody-lysozyme example: ready for prediction.")
        self.details.set("")
        self.errors.set("")
        self.read_structure()
        self._task_changed()
        self._review_state()
        action = "Choose a local model to predict." if self._needs_local_model() else "Select Predict to continue."
        self.status.set("1BVK is ready: antibody chains A,B and lysozyme chain C. " + action)

    def retry_layout(self):
        if not self._can_compute():
            return
        if self.document is None:
            self.errors.set("Open an analysis with a saved native interface first.")
            return
        path = self.document.path
        try:
            timeout = float(self.vars["layout_timeout"].get())
            if not math.isfinite(timeout) or timeout <= 0:
                raise ValueError("Layout timeout must be finite and positive")
        except ValueError as exc:
            self.errors.set(str(exc))
            return
        binary = self.vars["optcuts_bin"].get().strip() or "OptCuts_bin"
        self.status.set("Regenerating atlases from the saved interface; existing predictions are retained.")
        def action():
            from topoppi_affinity.workflow import layout_existing
            layout_existing(path, optcuts_bin=binary, layout_timeout=timeout, cancel_event=self.cancel,
                            log=lambda text: self.events.put(("log", text)),
                            event_callback=lambda event: self.events.put(("stage", event)))
            self.events.put(("result", path))
        self._launch(action)

    def _poll(self):
        self._poll_id = None
        if not self.window.winfo_exists():
            return
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "log":
                    self._append_activity(value)
                    self.status.set(value.splitlines()[-1] if value else "")
                elif kind == "stage":
                    if value.get("analysis_json") and value.get("stage") in {"prediction_ready", "layout_ready"} and not self.closing:
                        self.load(value["analysis_json"])
                    if value.get("message"):
                        self.status.set(value["message"])
                        self._append_activity(value["message"])
                    if value.get("total"):
                        self.progress.stop()
                        self.progress.configure(mode="determinate", maximum=value["total"], value=value.get("current", 0))
                    elif value.get("status") == "running":
                        self.progress.configure(mode="indeterminate")
                        self.progress.start()
                elif kind == "download":
                    message = value.get("message", "Preparing resources...")
                    if value.get("total") and value.get("current") is not None:
                        current, total = value["current"], value["total"]
                        message = f"{value.get('resource', 'Resource')}: {current / 1e6:.1f} / {total / 1e6:.1f} MB ({current / total:.0%})"
                    self.status.set(message)
                    if value.get("total"):
                        self.progress.stop()
                        self.progress.configure(mode="determinate", maximum=value["total"], value=value.get("current", 0))
                elif kind == "result" and not self.closing:
                    self.load(value)
                elif kind == "resources_ready":
                    self._task_changed()
                    self.status.set(value)
                elif kind in {"error", "cancelled"}:
                    self._append_activity(value)
                    self.status.set("Prediction stopped. Update the inputs or resources, then try again." if kind == "error" else value)
                    self.errors.set(value)
                    if kind == "error" and any(word in value for word in ("ModuleNotFoundError", "ImportError")):
                        self.errors.set(value + "\nOpen Advanced > Manage prediction dependencies to install the required libraries.")
                elif kind == "finished":
                    self._busy(False)
        except queue.Empty:
            pass
        if self.closing and self.close():
            return
        self._poll_id = self.window.after(100, self._poll)

    def open_result(self):
        path = filedialog.askopenfilename(parent=self.window, title="Open affinity analysis.json", filetypes=[("Analysis", "*.json")])
        if path:
            self.load(path)

    def load(self, path):
        from topoppi_affinity.workflow import load_analysis
        try:
            document = load_analysis(path)
            previous_selection = self.table.selection()
            self.document = document
            self.table.delete(*self.table.get_children())
            for row in document.predictions.itertuples():
                try:
                    number = float(row.predicted)
                    value = f"{number:.4g}" if math.isfinite(number) else "Unavailable"
                except (ValueError, TypeError):
                    value = "Unavailable"
                self.table.insert("", tk.END, iid=row.variant_id,
                                  values=(row.variant_id, row.mutations, value,
                                          {"base_only": "Base model only", "ok": "Complete", "error": "Failed", "pending": "Pending", "cancelled": "Cancelled"}.get(row.status, row.status)))
            report = document.report
            target = report.get("model", {}).get("target")
            units = {"pkd": "pKd (whole complex)", "ddg": "ddG kcal/mol (positive: weaker binding)"}.get(target, "Prediction target not yet determined")
            self.table.heading("prediction", text={"pkd": "pKd", "ddg": "ddG (kcal/mol)"}.get(target, "Prediction"))
            groups = report.get("interface", {}).get("groups", [[], []])
            self.partner_labels = [f"{name}: {','.join(group) or 'unavailable'}" for name, group in zip(("A", "B"), groups, strict=True)]
            self.partner_selector.configure(values=self.partner_labels)
            self.partner_label.set(self.partner_labels[int(self.side.get())])
            assignment = report.get("scope_assignment", {})
            scope = report.get("model", {}).get("scope", "Model scope unavailable")
            if assignment:
                scope += "\nSystem context: " + (assignment.get("user_basis") or "User-selected model scope")
            base_only = "status" in document.predictions and document.predictions.status.eq("base_only").any()
            base_note = "\nBASE MODEL ONLY: the required prior is unavailable; the displayed estimate uses the base component." if base_only else ""
            atlas_states = {"ok": "Ready", "not_requested": "Choose Generate 2D atlases to add maps",
                            "pending": "Waiting", "error": "Generation failed", "cancelled": "Cancelled"}
            layout_state = "Atlas: " + atlas_states.get(report['layout_status'], report['layout_status'])
            if report.get("layout_attempt_status"):
                attempt = report['layout_attempt_status']
                layout_state += " | Latest atlas attempt: " + atlas_states.get(attempt, attempt)
            source_path = report.get("input", {}).get("path")
            result_context = f"Displayed result: {Path(source_path).name}\n" if source_path else ""
            self.summary.set(f"{result_context}{units}\n{self.partner_labels[0]} | {self.partner_labels[1]}\n{scope}{base_note}\n{layout_state}")
            failures = [report[key] for key in ("interface_error", "prediction_error") if report.get(key)]
            if base_only and report.get("prior", {}).get("error"):
                failures.append(str(report["prior"]["error"]))
            notices = [f"Notice: {warning}" for warning in report.get("warnings", [])]
            self.errors.set("\n".join(dict.fromkeys(failures + notices)))
            completed = document.predictions.status.isin(["ok", "base_only"]).sum()
            review_action = ("Select a row to review its interface." if document.interface is not None else
                             "Select a row for status and error details; use Reuse inputs to prepare another run.")
            self.status.set(f"{completed} of {len(document.predictions)} results available. {review_action}")
            self._review_state()
            children = self.table.get_children()
            if children:
                self.table.selection_set(previous_selection[0] if previous_selection and previous_selection[0] in children else children[0])
            self.draw()
            self._append_activity(f"Opened analysis: {document.path}")
            if self.host is not None and hasattr(self.host, "receive_affinity_result"):
                self.host.receive_affinity_result(document)
        except (OSError, ValueError, KeyError) as exc:
            self.status.set(f"Cannot open result: {exc}")
            self.errors.set(f"Cannot open result: {exc}")

    def _choose_partner(self, _event=None):
        self.side.set(str(self.partner_labels.index(self.partner_label.get())))
        self.draw()

    def _table_selected(self, _event=None):
        if self.table.selection() != self._drawn_selection:
            self.draw()

    def draw(self):
        if self.document is None:
            return
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        selection = self.table.selection()
        self._drawn_selection = selection
        variant = selection[0] if selection else None
        side = int(self.side.get())
        self.partner_label.set(self.partner_labels[side])
        self.view_label.set(next(label for label, value in VIEWS.items() if value == self.mode.get()))
        self.descriptor_label.set(next(label for label, value in DESCRIPTORS.items() if value == self.descriptor.get()))
        detail = "Colours show native contact counts or areas (square angstroms). Predictions apply to the whole complex or variant."
        try:
            if variant:
                row, visible, hidden = self.document.selection(variant, side)
                total = "Joint variant total; highlighted residues are one group." if "," in row.mutations else "Effect of this specific substitution."
                if not row.mutations:
                    total = "Whole-complex affinity."
                coverage_names = {"absent_from_structure": "absent from the structure",
                                  "outside_retained_interface": "outside the retained interface",
                                  "footprint": "on the interface"}
                coverage = "; ".join(f"{site['chain']}:{site['position']}{site['insertion']} "
                    f"({coverage_names.get(site['coverage'], site['coverage']) if site['side'] in (side, None) else 'other partner'})"
                    for site in hidden)
                extra_labels = {"predicted_kd_M": "Predicted Kd (M)", "predicted_dg_at_reference_kcal_mol": "dG at reference temperature (kcal/mol)",
                                "anchored_mutant_kd_M": "Anchored mutant Kd (M)", "components_used": "Model components"}
                extras = [f"{label}: {row[name]}" for name, label in extra_labels.items() if name in row.index and row[name] != ""]
                sites = f" Substitution sites on this partner: {len(visible)}. {coverage}" if row.mutations else ""
                detail = (f"{total}{sites}\n" + "; ".join(extras) +
                          (f"\nPrediction: {row.error}" if row.error else "") + "\n" + detail)
            self.details.set(detail)
            patches, viz, style = self.document.view(side, self.mode.get(), variant,
                None if self.descriptor.get() == "none" else self.descriptor.get())
            if not patches:
                raise ValueError("No retained interface on this side; predictions remain in the table")
            figure = viz.plot_patches(patches, show=False, style_config=style)
            self._clear_plot()
            self.figure = figure
            self.canvas = FigureCanvasTkAgg(figure, master=self.plot)
            toolbar = NavigationToolbar2Tk(self.canvas, self.plot, pack_toolbar=False)
            toolbar.update()
            toolbar.pack(side=tk.BOTTOM, fill=tk.X)
            self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self.canvas.draw_idle()
        except (ValueError, KeyError, OSError) as exc:
            self._clear_plot()
            action = "Select 3D interface or Generate 2D atlases to continue." if self.mode.get() == "atlas" else "Predictions and residue coverage are available in the table."
            self.details.set(f"{detail}\nView: {exc}\n{action}")

    def export_predictions(self):
        if self.document is None:
            self.errors.set("Open or complete an analysis before exporting results.")
            return
        path = filedialog.asksaveasfilename(parent=self.window, title="Export predictions with units and status",
            initialfile="predictions.csv", defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if path:
            try:
                self.document.predictions.to_csv(path, index=False)
                if self.errors.get().startswith("Cannot export predictions:"):
                    self.errors.set(self._export_previous_error)
                self.status.set(f"Exported predictions, units, component status and errors to {path}")
            except OSError as exc:
                if not self.errors.get().startswith("Cannot export predictions:"):
                    self._export_previous_error = self.errors.get()
                self.errors.set(f"Cannot export predictions: {exc}")

    def open_folder(self):
        if self.document is None:
            self.errors.set("Open or complete an analysis first.")
            return
        path = str(self.document.path.parent)
        try:
            if sys.platform == "win32":
                os.startfile(path)
            else:
                subprocess.Popen(["open" if sys.platform == "darwin" else "xdg-open", path],
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except OSError as exc:
            self.errors.set(f"Cannot open the folder: {exc}. Results are saved at {path}")

    def reuse_inputs(self):
        if self.document is None:
            return
        report = self.document.report
        inputs = report.get("input", {})
        self.vars["structure"].set(inputs.get("path", ""))
        for key in ("chains_a", "chains_b", "source", "complex_id"):
            self.vars[key].set(inputs.get(key) or "")
        self.vars["temperature"].set(str(inputs.get("temperature", 298.15)))
        self.vars["wt_kd"].set(str(inputs["wt_kd"]) if inputs.get("wt_kd") is not None else "")
        self.vars["model"].set(report.get("model_path", ""))
        self.vars["output"].set(str(self.document.path.parent))
        for key in ("prior_csv", "prior_repository", "prior_checkpoint", "sequence_dir"):
            self.vars[key].set("")
        saved_prior = report.get("files", {}).get("prior")
        if saved_prior and (self.document.path.parent / saved_prior).is_file():
            self.vars["prior_csv"].set(str(self.document.path.parent / saved_prior))
        kind = report.get("model", {}).get("kind")
        from topoppi_affinity.resources import resource_status
        released = resource_status(cache_dir=self.vars["resource_cache"].get().strip() or None)["models"]
        recorded = report.get("model_path", "")
        uses_release = bool(recorded) and any(
            Path(recorded).expanduser() == Path(model["path"]).expanduser() for model in released.values())
        target = report.get("model", {}).get("target")
        self.task.set("Binding affinity (pKd)" if target == "pkd" or kind == "absolute" else "Mutation effects (ddG)")
        self.model_source.set("Prepared task model" if uses_release else "Local model file")
        context = report.get("scope_assignment", {}).get("context")
        self.mutation_context.set(context if context in {"known", "unseen"} else kind if kind in {"known", "unseen"} else "")
        self.vars["scope_note"].set(report.get("scope_assignment", {}).get("user_basis", ""))
        self.variants.delete("1.0", tk.END)
        rows = self.document.predictions
        if rows.mutations.ne("").any():
            self.variants.insert("1.0", "\n".join(f"{row.variant_id}\t{row.mutations}" for row in rows.itertuples()))
        self.retry.set(True)
        self._task_changed()
        if self.vars["structure"].get():
            self.read_structure()
        self.tabs.select(0)
        self.status.set("Inputs restored for a new run. Previous results will be kept in their original directory.")

    def show_in_main(self):
        if self.document is None or self.host is None or not hasattr(self.host, "receive_affinity_result"):
            return
        if self.host._busy:
            self.errors.set("The main window is computing an interface. Show this result when that run finishes.")
            return
        selection = self.table.selection()
        self.host.receive_affinity_result(self.document, variant_id=selection[0] if selection else None,
                                          side=int(self.side.get()), show=True)
        self.host.root.deiconify()
        self.host.root.lift()

    def _append_activity(self, text):
        line = f"[{datetime.now():%H:%M:%S}] {text}"
        self.log_history.append(line)
        self.log_history = self.log_history[-500:]
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, line + "\n")
        if int(self.log_text.index("end-1c").split(".")[0]) > 2000:
            self.log_text.delete("1.0", "500.0")
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)

    def diagnostics_text(self):
        data = {"python": sys.version, "executable": sys.executable, "task": self.task.get(),
                "mutation_context": self.mutation_context.get(), "model_source": MODEL_SOURCES[self.model_source.get()],
                "fields": {key: value.get() for key, value in self.vars.items()},
                "analysis": self.document.report if self.document is not None else None}
        return json.dumps(data, indent=2, default=str) + "\n\nActivity\n" + "\n".join(self.log_history)

    def copy_diagnostics(self):
        self.window.clipboard_clear()
        self.window.clipboard_append(self.diagnostics_text())
        self.status.set("Prediction setup and activity copied to the clipboard.")

    def save_diagnostics(self):
        path = filedialog.asksaveasfilename(parent=self.window, title="Save prediction diagnostics",
            initialfile="affinity-diagnostics.txt", defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if path:
            try:
                Path(path).write_text(self.diagnostics_text(), encoding="utf-8")
                self.status.set(f"Diagnostics saved to {path}")
            except OSError as exc:
                self.errors.set(f"Cannot save diagnostics: {exc}")

    def save_setup(self):
        """Save unfinished form text without loading a model or validating a run."""
        path = filedialog.asksaveasfilename(parent=self.window, title="Save prediction setup",
            initialfile="prediction.topoppi-setup.json", defaultextension=".topoppi-setup.json",
            filetypes=[("Prediction setup", "*.topoppi-setup.json"), ("JSON", "*.json")])
        if not path:
            return
        setup = self.export_setup()
        setup.pop("analysis", None)
        try:
            Path(path).write_text(json.dumps({"format": "topoppi-affinity-setup", "version": 1, "setup": setup}, indent=2) + "\n",
                                  encoding="utf-8")
        except OSError as exc:
            self.errors.set(f"Cannot save prediction setup: {exc}")
            return
        if self.errors.get().startswith("Cannot save prediction setup:"):
            self.errors.set("")
        self.status.set(f"Setup saved to {path}. Open setup restores these fields and unfinished variants in a later session.")

    def open_setup(self):
        """Opening a selected setup explicitly replaces the editable form."""
        if self._working or self.closing:
            self.errors.set("Wait for the active prediction to finish or cancel it before opening another setup.")
            return
        path = filedialog.askopenfilename(parent=self.window, title="Open prediction setup (replaces current setup)",
            filetypes=[("Prediction setup", "*.topoppi-setup.json"), ("JSON", "*.json")])
        if not path:
            return
        try:
            document = json.loads(Path(path).read_text(encoding="utf-8"))
            if not isinstance(document, dict) or document.get("format") != "topoppi-affinity-setup" or document.get("version") != 1:
                raise ValueError("Choose a prediction setup created with Save setup. To review an analysis, use Open saved result.")
            setup = document.get("setup")
            if (not isinstance(setup, dict) or not isinstance(setup.get("fields"), dict)
                    or not all(isinstance(value, str) for value in setup["fields"].values())
                    or not isinstance(setup.get("variants"), str)
                    or setup.get("prediction_task") not in ("absolute", "mutation")
                    or setup.get("mutation_context") not in ("", "known", "unseen")
                    or setup.get("model_source") not in ("prepared", "local")
                    or not isinstance(setup.get("layout"), bool) or not isinstance(setup.get("retry"), bool)):
                raise ValueError("The file has invalid setup fields. Choose a file created with Save setup.")
            # Setup files contain form values, not a dependency on an old result.
            form = {key: setup[key] for key in ("fields", "prediction_task", "mutation_context", "model_source",
                                               "variants", "layout", "retry")}
            if self.errors.get().startswith("Cannot open prediction setup:"):
                self.errors.set("")
            self.restore_setup(form)
            self.status.set(f"Setup opened from {path}. Review source paths and unfinished variants before predicting.")
        except (OSError, ValueError, RuntimeError) as exc:
            self.errors.set(f"Cannot open prediction setup: {exc}")

    def binding_setup(self):
        """Share structure identity; model protocols and variant drafts stay in this panel."""
        return {name: self.vars[name].get().strip() for name in ("structure", "chains_a", "chains_b")}

    def export_setup(self):
        """Preserve the actual edited form across an explicit dependency restart."""
        return {"fields": {name: variable.get() for name, variable in self.vars.items()},
                "task": self.task.get(), "prediction_task": TASKS[self.task.get()],
                "mutation_context": self.mutation_context.get(), "model_source": MODEL_SOURCES[self.model_source.get()],
                "variants": self.variants.get("1.0", "end-1c"),
                "layout": self.layout.get(), "retry": self.retry.get(),
                "analysis": str(self.document.path) if self.document is not None else None}

    def _check_shared_source(self, shared, *, initial=False):
        from topoppi.interface import chain_group

        if initial or not self.variants.get("1.0", "end-1c").strip():
            return
        current = self.binding_setup()
        different_structure = (current["structure"] and
            Path(current["structure"]).expanduser().resolve() != Path(shared["structure"]).expanduser().resolve())
        different_partners = (current["structure"] and any(
            set(chain_group(current[key])) != set(chain_group(shared[key])) for key in ("chains_a", "chains_b")))
        if different_structure or different_partners:
            raise ValueError("The mutation draft uses a different structure or partner groups. Your edits were kept. "
                             "Open Mutation effects, save any variant draft you need, and explicitly switch its structure and partners before adding these residues.")

    def _append_variant_setup(self, setup, *, initial=False):
        """Validate a residue handoff before changing any edited form field."""
        import pandas as pd

        from topoppi_affinity.workflow import inspect_structure, variants_table

        incoming = variants_table(variants=pd.DataFrame(setup["append_variants"]))
        if incoming.empty:
            raise ValueError("Select at least one residue and define a substitution before creating variants.")
        shared = {key: str(setup.get(key, "")).strip() for key in ("structure", "chains_a", "chains_b")}
        if not all(shared.values()):
            raise ValueError("Variant creation needs the source structure and both partner groups.")
        self._check_shared_source(shared, initial=initial)
        records = inspect_structure(shared["structure"])
        previous_records = self.structure_records
        try:
            self.structure_records = records
            for index in incoming.index:
                self._check_structure_inputs({**shared, "task": "unseen", "variants": incoming.loc[[index]]})
        finally:
            self.structure_records = previous_records
        existing_text = self.variants.get("1.0", "end-1c")
        try:
            existing = variant_lines(existing_text)
        except ValueError as exc:
            raise ValueError(f"Correct the existing variant rows before appending residues: {exc}. Your edits were kept.") from exc
        used = set(existing.variant_id) if existing is not None else set()
        rows = []
        for row in incoming.itertuples():
            identifier, suffix = row.variant_id, 2
            while identifier in used:
                identifier = f"{row.variant_id}_{suffix}"
                suffix += 1
            used.add(identifier)
            rows.append(f"{identifier}\t{row.mutations}")
        for key, value in shared.items():
            self.vars[key].set(value)
        self.task.set("Mutation effects (ddG)")
        separator = "\n" if existing_text and not existing_text.endswith("\n") else ""
        self.variants.insert(tk.END, separator + "\n".join(rows))
        self.read_structure(records=records)
        self.tabs.select(0)
        self.errors.set("")
        self.status.set(f"Added {len(rows)} variant{'s' if len(rows) != 1 else ''} from the selected residues. "
                        "Existing rows and model settings are retained. Review model applicability before prediction.")
        self.variants.see(tk.END)

    def restore_setup(self, setup, *, initial=False):
        if self._working or self.closing or (self.worker is not None and self.worker.is_alive()):
            raise RuntimeError("A binding prediction is active. Wait for it to finish or cancel it before changing tasks or adding variants.")
        if "append_variants" in setup:
            self._append_variant_setup(setup, initial=initial)
            return
        if any(key in setup for key in ("structure", "chains_a", "chains_b")):
            shared = {key: str(setup.get(key, value)).strip() for key, value in self.binding_setup().items()}
            self._check_shared_source(shared, initial=initial)
        for name, value in {**setup.get("fields", {}), **{key: setup[key] for key in
                            ("structure", "chains_a", "chains_b") if key in setup}}.items():
            if name in self.vars:
                self.vars[name].set(value)
        # Existing restart files used deployment-task names as the task selector.
        task = setup.get("prediction_task", setup.get("task", TASKS[self.task.get()]))
        legacy_scopes = {"Mutation effects: known system": "known", "Mutation effects: unseen system": "unseen"}
        task = legacy_scopes.get(task, TASKS.get(task, task))
        if task in {"known", "unseen"}:
            self.mutation_context.set(task)
            task = "mutation"
        if task == "Custom model":
            self.model_source.set("Local model file")
        else:
            self.task.set(next(label for label, key in TASKS.items() if key == task))
        if "mutation_context" in setup:
            self.mutation_context.set(setup["mutation_context"])
        if "model_source" in setup:
            self.model_source.set(next(label for label, key in MODEL_SOURCES.items() if key == setup["model_source"]))
        if "variants" in setup:
            self.variants.delete("1.0", tk.END)
            self.variants.insert("1.0", setup["variants"])
        self.layout.set(setup.get("layout", self.layout.get()))
        self.retry.set(setup.get("retry", self.retry.get()))
        self._task_changed()
        if setup.get("analysis"):
            self.load(setup["analysis"])
        else:
            if self.vars["structure"].get():
                self.read_structure()
            self.tabs.select(0)
            self.status.set("Prediction setup restored. Review the inputs, then select Predict.")

    def _clear_plot(self):
        import matplotlib.pyplot as plt

        if self.canvas is not None and self.canvas._idle_draw_id:
            self.canvas.get_tk_widget().after_cancel(self.canvas._idle_draw_id)
            self.canvas._idle_draw_id = None
        for child in self.plot.winfo_children():
            child.destroy()
        if self.figure is not None:
            plt.close(self.figure)
        self.figure, self.canvas = None, None

    def close(self):
        """Cooperative shutdown; Tk stays responsive until the worker stops."""
        if not self.window.winfo_exists():
            return True
        self.closing = True
        self.cancel.set()
        installer = getattr(self, "_affinity_installer", None)
        if installer is not None and installer.window.winfo_exists() and not installer.close():
            self.status.set("Closing after the dependency installer stops...")
            return False
        if self.worker is not None and self.worker.is_alive():
            self.status.set("Closing after the active calculation stops...")
            return False
        if self.host is not None:
            self.host._affinity_draft = self.export_setup()
        if self._poll_id is not None:
            self.window.after_cancel(self._poll_id)
            self._poll_id = None
        if self.figure is not None or self.canvas is not None:
            self._clear_plot()
        self.window.destroy()
        return True


def open_panel(host):
    previous = getattr(host, "_affinity_panel", None)
    if previous is not None and previous.window.winfo_exists():
        previous._refresh_main_input_context()
        previous.window.lift()
        return previous
    host._affinity_panel = AffinityPanel(host.root, host=host)
    if draft := getattr(host, "_affinity_draft", None):
        host._affinity_panel.restore_setup(draft)
        host._affinity_panel.status.set("Your unfinished prediction setup was restored from this session. Review it before predicting.")
    return host._affinity_panel


def main(analysis=None):
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise RuntimeError("Cannot open a Tk graphical session. Start the GUI on a desktop, or use "
                           "topoppi-affinity analyze / inspect / render from the command line.") from exc
    root.withdraw()
    panel = AffinityPanel(root, analysis=analysis)
    def finish():
        if panel.window.winfo_exists():
            root.after(100, finish)
        else:
            root.destroy()
    root.after(100, finish)
    root.mainloop()
