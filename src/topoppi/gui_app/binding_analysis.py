"""Shared task navigation and residue-to-variant handoff for optional analyses."""
from __future__ import annotations

import re
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter import font as tkfont

from Bio.Data.PDBData import protein_letters_3to1_extended

BINDING_TASKS = {
    "absolute": ("Binding affinity", "Predict absolute affinity (pKd / Kd).", "affinity"),
    "mutation": ("Mutation effects", "Predict binding changes for explicit variants (ddG).", "affinity"),
    "hotspot": ("Rank hotspots", "Rank residues within one selected partner.", "hotspot"),
}
AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")


def available_variant_structure(path):
    path = Path(path).expanduser()
    return path.suffix.lower() in {".pdb", ".cif", ".mmcif", ".ent"} and path.is_file()


def validate_variant_structure(path, setup, residues):
    """Check a user-selected replacement against the author identities being transferred."""
    from Bio.PDB.PDBExceptions import PDBConstructionException

    from topoppi.interface import chain_group
    from topoppi.structure import load_structure

    if not available_variant_structure(path):
        raise ValueError("Choose an existing PDB or mmCIF structure file.")
    try:
        structure = load_structure(path)
    except (OSError, ValueError, KeyError, PDBConstructionException) as exc:
        raise ValueError(f"Could not read this structure. Choose a readable PDB/mmCIF file. Details: {exc}") from exc
    groups = [set(chain_group(setup[key])) for key in ("chains_a", "chains_b")]
    if groups[0] & groups[1]:
        raise ValueError("Assign each chain to one binding partner; the partner groups overlap.")
    missing = (groups[0] | groups[1]) - set(structure.chain_ids)
    if missing:
        raise ValueError("This file is missing the recorded partner chains: " + ", ".join(sorted(missing)) + ". Choose the matching source structure.")
    identities = {}
    for chain in structure.chains:
        for atom in chain.atoms:
            residue = atom.get_parent()
            _, number, insertion = residue.get_id()
            identities[(chain.chain_id, number, insertion.strip())] = protein_letters_3to1_extended.get(residue.get_resname().upper(), "X")
    for row in residues.values():
        key = (row["chain"], row["position"], row["insertion"])
        if row["chain"] not in groups[0] | groups[1] or identities.get(key) != row["wt"]:
            label = f"{row['chain']}:{row['wt']}{row['position']}{row['insertion']}"
            raise ValueError(f"This file does not match selected residue {label}. Choose the matching source structure.")
    return str(Path(path).expanduser().resolve())


def normalize_residues(residues):
    """Preserve author numbering; never infer identities from display labels."""
    result = {}
    for row in residues:
        chain, insertion, wt = (str(row.get(name, "")).strip() for name in ("chain", "insertion", "wt"))
        position = str(row.get("position", "")).strip()
        if (not re.fullmatch(r"[^:\s,]+", chain) or not re.fullmatch(r"-?\d+", position)
                or not re.fullmatch(r"[A-Z]?", insertion) or wt not in AMINO_ACIDS):
            raise ValueError("Selected residues need a chain, author residue number, insertion code and standard amino acid.")
        key = f"{chain}:{int(position)}{insertion}"
        value = {"chain": chain, "position": int(position), "insertion": insertion, "wt": wt}
        if key in result and result[key] != value:
            raise ValueError(f"Conflicting identities for residue {key}.")
        result[key] = value
    if not result:
        raise ValueError("Select at least one residue first.")
    return result


def variant_rows(residues, targets, mode):
    mutations = []
    for key, row in residues.items():
        target = targets[key]
        if target not in AMINO_ACIDS:
            raise ValueError(f"Choose a target amino acid for {key}.")
        if target == row["wt"]:
            raise ValueError(f"{key} is already {target}. Choose a different amino acid or remove this residue.")
        mutations.append(f"{row['chain']}:{row['wt']}{row['position']}{row['insertion']}{target}")
    if mode == "joint":
        return [{"variant_id": "joint_variant", "mutations": ",".join(mutations)}]
    return [{"variant_id": value.replace(":", "_"), "mutations": value} for value in mutations]


class VariantComposer:
    def __init__(self, host, residues, setup):
        self.host, self.setup = host, dict(setup)
        self.residues = normalize_residues(residues)
        self.window = tk.Toplevel(host.root)
        self.window.title("Create variants - Binding analysis")
        self.window.geometry("700x600")
        self.window.minsize(620, 480)
        body = ttk.Frame(self.window, padding=14)
        body.pack(fill=tk.BOTH, expand=True)
        body.columnconfigure(0, weight=1)
        body.rowconfigure((3, 6), weight=1, minsize=45)
        ttk.Label(body, text="Create variants from selected residues", font=("TkDefaultFont", 12, "bold")).grid(row=0, sticky=tk.W)
        source = ttk.Frame(body)
        source.grid(row=1, sticky=tk.EW, pady=(5, 10))
        source.columnconfigure(0, weight=1)
        self.source_path = tk.StringVar(self.window, str(setup.get("structure", "")))
        ttk.Entry(source, textvariable=self.source_path, state="readonly", width=45).grid(row=0, column=0, sticky=tk.EW)
        ttk.Button(source, text="Locate structure...", command=self.locate_structure).grid(row=0, column=1, padx=(6, 0))
        ttk.Label(body, text=f"Partners: {setup['chains_a']} / {setup['chains_b']}. Select a residue and its target amino acid (default: A).",
                  wraplength=560).grid(row=2, sticky=tk.W)
        table_frame = ttk.Frame(body)
        table_frame.grid(row=3, sticky=tk.NSEW, pady=6)
        self.table = ttk.Treeview(table_frame, columns=("residue", "target"), show="headings", height=5, selectmode="browse")
        self.table.heading("residue", text="Selected residue (author numbering)")
        self.table.heading("target", text="Change to")
        self.table.column("residue", width=330)
        self.table.column("target", width=130)
        table_scroll = ttk.Scrollbar(table_frame, command=self.table.yview)
        table_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.table.configure(yscrollcommand=table_scroll.set)
        self.table.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.targets = {key: tk.StringVar(self.window, "A") for key in self.residues}
        for key, row in self.residues.items():
            self.table.insert("", "end", iid=key, values=(f"{row['chain']}:{row['wt']}{row['position']}{row['insertion']}", "A"))
            self.targets[key].trace_add("write", lambda *_: self.refresh_preview())
        edit = ttk.Frame(body)
        edit.grid(row=4, sticky=tk.EW)
        ttk.Label(edit, text="Target amino acid:").pack(side=tk.LEFT)
        self.target = tk.StringVar(self.window, "A")
        target = ttk.Combobox(edit, textvariable=self.target, values=AMINO_ACIDS, state="readonly", width=5)
        target.pack(side=tk.LEFT, padx=8)
        target.bind("<<ComboboxSelected>>", self._edit_target)
        ttk.Button(edit, text="Remove selected", command=self._remove).pack(side=tk.LEFT)
        self.table.bind("<<TreeviewSelect>>", self._select)
        self.mode = tk.StringVar(self.window, "single")
        modes = ttk.LabelFrame(body, text="How should these changes be grouped?", padding=8)
        modes.grid(row=5, sticky=tk.EW, pady=10)
        ttk.Radiobutton(modes, text="Separate single variants: one prediction per residue", variable=self.mode,
                        value="single", command=self.refresh_preview).pack(anchor=tk.W)
        ttk.Radiobutton(modes, text="One joint variant: all changes in the same protein complex", variable=self.mode,
                        value="joint", command=self.refresh_preview).pack(anchor=tk.W)
        preview_frame = ttk.Frame(body)
        preview_frame.grid(row=6, sticky=tk.NSEW)
        preview_frame.columnconfigure(0, weight=1)
        preview_frame.rowconfigure(0, weight=1)
        self.preview = ttk.Treeview(preview_frame, columns=("id", "mutations"), show="headings", height=4)
        self.preview.heading("id", text="Variant")
        self.preview.heading("mutations", text="Exact changes")
        self.preview.column("id", width=150, stretch=False)
        self.preview.column("mutations", width=400, stretch=False)
        preview_scroll = ttk.Scrollbar(preview_frame, command=self.preview.yview)
        preview_scroll.grid(row=0, column=1, sticky=tk.NS)
        preview_horizontal = ttk.Scrollbar(preview_frame, orient=tk.HORIZONTAL, command=self.preview.xview)
        preview_horizontal.grid(row=1, column=0, sticky=tk.EW)
        self.preview.configure(yscrollcommand=preview_scroll.set, xscrollcommand=preview_horizontal.set)
        self.preview.grid(row=0, column=0, sticky=tk.NSEW)
        self.preview.bind("<Control-c>", self.copy_variants)
        self.preview.bind("<Command-c>", self.copy_variants)
        self.feedback = tk.StringVar(self.window)
        ttk.Label(body, textvariable=self.feedback, wraplength=560, foreground="#52616B").grid(row=7, sticky=tk.W, pady=8)
        actions = ttk.Frame(body)
        actions.grid(row=8, sticky=tk.EW)
        self.copy_button = ttk.Button(actions, text="Copy variants", command=self.copy_variants)
        self.copy_button.pack(side=tk.LEFT)
        self.add_button = ttk.Button(actions, text="Add to mutation analysis", command=self.submit, style="Primary.TButton")
        self.add_button.pack(side=tk.RIGHT)
        ttk.Button(actions, text="Cancel", command=self.window.destroy).pack(side=tk.RIGHT, padx=8)
        self.table.selection_set(next(iter(self.residues)))
        self.refresh_preview()

    def locate_structure(self):
        path = filedialog.askopenfilename(parent=self.window, title="Locate the source structure for these residues",
                                         filetypes=[("Protein structures", "*.pdb *.cif *.mmcif *.ent")])
        if not path:
            return
        try:
            source = validate_variant_structure(path, self.setup, self.residues)
        except (ValueError, OSError) as exc:
            self.feedback.set(str(exc))
            return
        # This is an explicit input for the new mutation task, not a rewrite of
        # the displayed atlas or its recorded provenance.
        self.setup["structure"] = source
        self.source_path.set(source)
        self.refresh_preview()

    def _select(self, _event=None):
        selected = self.table.selection()
        if selected:
            self.target.set(self.targets[selected[0]].get())

    def _edit_target(self, _event=None):
        selected = self.table.selection()
        if selected:
            self.targets[selected[0]].set(self.target.get())

    def _remove(self):
        for key in self.table.selection():
            self.table.delete(key)
            del self.residues[key]
            del self.targets[key]
        if self.residues:
            self.table.selection_set(next(iter(self.residues)))
        self.refresh_preview()

    def refresh_preview(self):
        self.preview.delete(*self.preview.get_children())
        self.copy_button.configure(text="Copy variants", state=tk.DISABLED)
        for key, variable in self.targets.items():
            self.table.set(key, "target", variable.get())
        try:
            if not self.residues:
                raise ValueError("No residues remain. Cancel and select residues to create variants.")
            rows = variant_rows(self.residues, {key: value.get() for key, value in self.targets.items()}, self.mode.get())
        except ValueError as exc:
            self.feedback.set(str(exc))
            self.add_button.configure(state=tk.DISABLED)
            return []
        for row in rows:
            self.preview.insert("", "end", values=(row["variant_id"], row["mutations"]))
        font_name = ttk.Style(self.window).lookup("Treeview", "font") or "TkDefaultFont"
        font = tkfont.Font(root=self.window, font=font_name)
        for column, field, minimum in (("id", "variant_id", 150), ("mutations", "mutations", 400)):
            self.preview.column(column, width=max(minimum, max(font.measure(row[field]) for row in rows) + 20))
        self.copy_button.configure(state=tk.NORMAL)
        if not available_variant_structure(self.setup.get("structure", "")):
            self.feedback.set("The recorded source file is unavailable. Choose Locate structure to supply the matching PDB/mmCIF. Your selected residues and changes are kept.")
            self.add_button.configure(state=tk.DISABLED)
            return []
        self.feedback.set(f"{len(rows)} {'variant' if len(rows) == 1 else 'variants'} will be added. Review them before running ddG prediction.")
        self.add_button.configure(state=tk.NORMAL)
        return rows

    def copy_variants(self, _event=None):
        keys = self.preview.selection() or self.preview.get_children()
        if keys:
            lines = ["variant_id\tmutations"]
            lines.extend("\t".join(map(str, self.preview.item(key, "values"))) for key in keys)
            self.window.clipboard_clear()
            self.window.clipboard_append("\n".join(lines))
            self.copy_button.configure(text=f"Copied {len(keys)} {'variant' if len(keys) == 1 else 'variants'}")
        return "break"

    def submit(self):
        rows = self.refresh_preview()
        if not rows:
            return
        try:
            self.host.open_binding_task("mutation", setup={**self.setup, "append_variants": rows})
        except (ValueError, RuntimeError, OSError) as exc:
            self.feedback.set(str(exc))
            return
        self.window.destroy()


class BindingAnalysisMixin:
    def _init_binding_analysis(self):
        if hasattr(self, "_binding_tasks"):
            return
        self._binding_tasks = {}
        self._binding_selection = {}
        self._binding_selection_artists = []
        self.var_select_binding_residues = tk.BooleanVar(self.root, False)
        self.binding_selection_label = tk.StringVar(self.root, "0 selected")
        self.btn_binding_analysis = ttk.Button(self.plugin_actions_frame, text="Binding analysis...", command=self.show_binding_analysis)
        self.btn_binding_analysis.pack(side=tk.LEFT, padx=(8, 0))
        self.analysis_menu.add_separator()
        self.analysis_menu.add_command(label="Binding analysis...", command=self.show_binding_analysis)
        for task, (title, _description, _plugin) in BINDING_TASKS.items():
            self.analysis_menu.add_command(label=title + "...", command=lambda task=task: self._navigate_binding_task(task))

    def register_binding_task(self, task, *, title, description, plugin, open_task):
        self._binding_tasks[task] = {"title": title, "description": description, "plugin": plugin, "open": open_task}

    def binding_setup(self, displayed=False):
        if displayed:
            run = self._successful_single_run
            if not run:
                return {}
            params = run["params"]
            geometry = run.get("geometry")
            interface = geometry.metadata() if geometry is not None else run.get("manifest", {}).get("interface", {})
            groups = interface.get("groups", [])
            if len(groups) == 2:
                return {"structure": interface.get("source", {}).get("path") or params.get("path", ""),
                        "chains_a": ",".join(groups[0]), "chains_b": ",".join(groups[1])}
            return {"structure": params.get("path", ""), "chains_a": params.get("chain_a", run["viz"].chain_a_id),
                    "chains_b": params.get("chain_b", run["viz"].chain_b_id)}
        return {"structure": self.var_input_path.get(), "chains_a": self.var_chain_a.get(), "chains_b": self.var_chain_b.get()}

    def open_binding_task(self, task, setup=None):
        if self._busy or self._closing:
            raise RuntimeError("Wait for the current analysis to finish before opening another task.")
        registered = self._binding_tasks.get(task)
        if registered is None:
            title, _description, plugin = BINDING_TASKS[task]
            raise RuntimeError(f"{title} needs the TopoPPI {plugin.title()} plugin. Open Plugins to install it, then refresh.")
        values = dict(self.binding_setup() if setup is None else setup)
        if task in {"absolute", "mutation"}:
            values["prediction_task"] = task
        return registered["open"](setup=values)

    def _navigate_binding_task(self, task, setup=None, parent=None):
        try:
            return self.open_binding_task(task, setup)
        except (ValueError, RuntimeError, OSError) as exc:
            messagebox.showerror("Binding analysis", str(exc), parent=parent or self.root)

    def add_binding_navigation(self, parent, current_task, setup_provider):
        frame = ttk.Frame(parent, padding=(10, 6))
        frame.pack(fill=tk.X)
        ttk.Label(frame, text="Tasks for current inputs").pack(side=tk.LEFT, padx=(0, 10))
        for task, (title, _description, _plugin) in BINDING_TASKS.items():
            def navigate(task=task):
                active = current_task() if callable(current_task) else current_task
                if task != active:
                    self._navigate_binding_task(task, setup_provider(), parent)
            ttk.Button(frame, text=title, command=navigate).pack(side=tk.LEFT, padx=3)
        return frame

    def show_binding_analysis(self):
        existing = getattr(self, "_binding_window", None)
        if existing is not None and existing.winfo_exists():
            self._binding_launcher_refresh()
            existing.deiconify()
            existing.lift()
            return
        window = self._binding_window = tk.Toplevel(self.root)
        window.title("Binding analysis")
        window.resizable(False, False)
        body = ttk.Frame(window, padding=18)
        body.pack(fill=tk.BOTH, expand=True)
        ttk.Label(body, text="What would you like to learn?", font=("TkDefaultFont", 13, "bold")).pack(anchor=tk.W)
        source = tk.StringVar(window, "inputs")
        source_options = ttk.Frame(body)
        source_options.pack(fill=tk.X, pady=10)
        summary = tk.StringVar(window)
        task_descriptions = {}
        def update_source():
            displayed_option.configure(state=tk.NORMAL if self._successful_single_run else tk.DISABLED)
            setup = self.binding_setup(source.get() == "displayed")
            summary.set(f"{setup.get('structure') or 'Choose a structure in the task'}\nPartners: {setup.get('chains_a') or '-'} / {setup.get('chains_b') or '-'}")
            for task, label in task_descriptions.items():
                _title, description, plugin = BINDING_TASKS[task]
                suffix = "" if task in self._binding_tasks else f"\nInstall the {plugin.title()} plugin to use this task."
                label.configure(text=description + suffix)
        ttk.Radiobutton(source_options, text="Current input fields", variable=source, value="inputs", command=update_source).pack(side=tk.LEFT)
        displayed_option = ttk.Radiobutton(source_options, text="Displayed result", variable=source, value="displayed", command=update_source)
        displayed_option.pack(side=tk.LEFT, padx=12)
        ttk.Label(body, textvariable=summary, wraplength=580, foreground="#52616B").pack(anchor=tk.W, pady=(0, 12))
        update_source()
        self._binding_launcher_refresh = update_source
        window.bind("<FocusIn>", lambda _event: update_source())
        for task, (title, _description, _plugin) in BINDING_TASKS.items():
            row = ttk.Frame(body)
            row.pack(fill=tk.X, pady=7)
            ttk.Button(row, text=title, width=23, command=lambda task=task: self._navigate_binding_task(
                task, self.binding_setup(source.get() == "displayed"), window)).pack(side=tk.LEFT, padx=(0, 12))
            task_descriptions[task] = ttk.Label(row, wraplength=360)
            task_descriptions[task].pack(side=tk.LEFT)
        update_source()
        ttk.Separator(body).pack(fill=tk.X, pady=12)
        ttk.Label(body, text="Select residues in a hotspot table or interface map, then choose Create variants.\nEach task keeps its own output units and model requirements.", wraplength=580).pack(anchor=tk.W)
        ttk.Button(body, text="Plugins...", command=self.show_plugins).pack(anchor=tk.E, pady=(12, 0))

    def create_variants_from_residues(self, residues, setup=None):
        if self._busy or self._closing:
            raise RuntimeError("Wait for the current analysis to finish before creating variants.")
        setup = dict(self.binding_setup(displayed=True) if setup is None else setup)
        if not setup.get("chains_a") or not setup.get("chains_b"):
            raise ValueError("The displayed result does not record both partner groups. Open its original structure and select both partners first.")
        composer = VariantComposer(self, residues, setup)
        return composer

    def _toggle_binding_residue(self, key):
        run = self._successful_single_run
        if not run or self._busy:
            return
        if key in self._binding_selection:
            del self._binding_selection[key]
        else:
            viz = run["viz"]
            metadata = {**viz.residue_metadata_B, **viz.residue_metadata_A}.get(key)
            if metadata is None:
                self.log("This residue has no recorded identity and cannot be used to create variants.")
                return
            residue_id = metadata["residue_id"]
            row = {"chain": metadata["chain"], "position": residue_id[1], "insertion": residue_id[2].strip(),
                   "wt": protein_letters_3to1_extended.get(metadata["residue_name"].upper(), "X")}
            try:
                normalize_residues([row])
            except ValueError as exc:
                self.log(str(exc))
                return
            self._binding_selection[key] = row
        self._show_binding_selection()

    def _show_binding_selection(self):
        changed = bool(self._binding_selection_artists)
        for artist in self._binding_selection_artists:
            artist.remove()
        self._binding_selection_artists = []
        self.binding_selection_label.set(f"{len(self._binding_selection)} selected")
        button = getattr(self, "btn_create_selected_variants", None)
        if button is not None and button.winfo_exists():
            button.configure(state=tk.NORMAL if self._binding_selection and not self._busy else tk.DISABLED)
        run = self._successful_single_run
        if run and self.current_fig:
            for objects in run["viz"].artist_map.values():
                if objects.get("residue_key") not in self._binding_selection:
                    continue
                anchor = self._label_anchor(objects)
                if anchor is not None:
                    reference = next((objects.get(name) for name in ("collection", "scatter", "text")
                                      if objects.get(name) is not None), None)
                    ax = reference.axes if reference is not None else self.current_fig.axes[0]
                    point = ax.scatter(*anchor, s=140, facecolors="none", edgecolors="#E58A18", linewidths=2.2, zorder=20,
                                       **({"depthshade": False} if self._surface_view() else {}))
                    self._binding_selection_artists.append(point)
            if changed or self._binding_selection_artists:
                self.current_canvas.draw_idle()

    def _clear_binding_selection(self):
        self._binding_selection.clear()
        self._show_binding_selection()

    def _create_selected_variants(self):
        try:
            return self.create_variants_from_residues(list(self._binding_selection.values()))
        except (ValueError, RuntimeError) as exc:
            messagebox.showerror("Create variants", str(exc), parent=self.root)
