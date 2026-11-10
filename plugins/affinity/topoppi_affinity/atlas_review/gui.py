"""Linked structure, contact graph and atlas review window. Graph and atlas consume the same native context."""

from __future__ import annotations

import json
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np

from .session import VIEW_CONDITIONS, VIEW_LABELS, ReviewSession


class ReviewWindow:
    def __init__(self, root, session, *, condition_locked=False, event_callback=None):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure

        self.session = session
        self.event_callback = event_callback or (lambda event: None)
        self.condition_locked = condition_locked
        self.window = tk.Toplevel(root)
        self.window.title("TopoPPI Affinity — review interface contacts")
        self.window.geometry(
            f"{min(1280, self.window.winfo_screenwidth() - 60)}x{min(850, self.window.winfo_screenheight() - 80)}"
        )
        controls = ttk.Frame(self.window)
        controls.pack(fill="x")
        self.query = tk.StringVar()
        ttk.Entry(controls, textvariable=self.query, width=24).pack(side="left")
        ttk.Button(controls, text="Find residues", command=self.search).pack(side="left")
        ttk.Button(controls, text="Clear selection", command=lambda: self.choose([])).pack(side="left")
        ttk.Button(controls, text="Save review...", command=self.save).pack(side="left")
        ttk.Button(controls, text="Open review...", command=self.open).pack(side="left")
        ttk.Button(controls, text="Export contacts...", command=self.export_contacts).pack(side="left")
        view_controls = ttk.Frame(self.window)
        view_controls.pack(fill="x")
        self.expanded = tk.BooleanVar(value=False)
        ttk.Checkbutton(view_controls, text="Enlarge 3D view", variable=self.expanded, command=self.draw).pack(side="left")
        ttk.Label(view_controls, text="View").pack(side="left", padx=(8, 4))
        self.condition = tk.StringVar(value=VIEW_LABELS[session.condition])
        box = ttk.Combobox(
            view_controls, textvariable=self.condition, values=tuple(VIEW_LABELS.values()), width=16, state="readonly"
        )
        box.pack(side="left")
        box.configure(state="disabled" if condition_locked else "readonly")
        box.bind("<<ComboboxSelected>>", lambda e: self.draw())
        self.status = tk.StringVar()
        ttk.Label(self.window, textvariable=self.status, wraplength=1200).pack(fill="x")
        self.variant = tk.StringVar(value=session.variant_id or "")
        ttk.Label(view_controls, text="Prediction").pack(side="left", padx=(8, 4))
        variants = ttk.Combobox(
            view_controls,
            textvariable=self.variant,
            values=list(session.document.predictions.variant_id.astype(str)),
            state="readonly",
            width=24,
        )
        variants.pack(side="left")
        variants.bind("<<ComboboxSelected>>", self.variant_selected)
        self.output = tk.StringVar(
            value="Prediction applies to the complete specified substitution set; surface colors show input context."
        )
        ttk.Label(self.window, textvariable=self.output, wraplength=1200).pack(fill="x")
        split = ttk.Panedwindow(self.window, orient="horizontal")
        split.pack(fill="both", expand=True)
        tableframe = ttk.Frame(split)
        plotframe = ttk.Frame(split)
        split.add(tableframe, weight=0)
        split.add(plotframe, weight=1)
        self.table = ttk.Treeview(tableframe, columns=("side", "area", "degree"), selectmode="extended", height=24)
        self.table.heading("#0", text="Residue (author ID)")
        self.table.column("#0", width=150)
        for name, label in [("side", "Partner"), ("area", "Area Å²"), ("degree", "Contacts")]:
            self.table.heading(name, text=label)
            self.table.column(name, width=65)
        self.table.tag_configure("neighbor", foreground="#228833")
        self.table.pack(fill="both", expand=True)
        self.syncing = False
        residues = session.document.residues
        area = (
            dict(zip(residues.residue_key, residues.footprint_area_A2, strict=True))
            if "footprint_area_A2" in residues
            else {}
        )
        for key in sorted(session.keys):
            self.table.insert(
                "",
                "end",
                iid=key,
                text=key,
                values=(session.side_by_key[key], f"{area.get(key, 0):.2f}", len(session.contacts.get(key, {}))),
            )
        self.table.bind("<<TreeviewSelect>>", self.table_selected)
        self.context_text = tk.Text(tableframe, height=9, width=44, wrap="word")
        self.context_text.pack(fill="x")
        self.measurements = []
        measurement_file = session.document.report.get("files", {}).get("measurements")
        if measurement_file:
            import pandas as pd

            self.measurements = pd.read_csv(
                session.document.path.parent / measurement_file, keep_default_na=False
            ).to_dict("records")
        ttk.Button(view_controls, text="Experimental measurements", command=self.show_measurements).pack(side="left")
        self.figure = Figure(figsize=(10, 7))
        self.canvas = FigureCanvasTkAgg(self.figure, master=plotframe)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        NavigationToolbar2Tk(self.canvas, plotframe)
        self.canvas.mpl_connect("pick_event", self.pick)
        self.canvas.mpl_connect("button_press_event", self.pick_surface)
        self.pick_maps = {}
        self.atlas_cache = {}
        self.draw()
        self.choose(session.selected)
        if session.variant_id:
            kept = set(session.selected)
            self.variant_selected(None)
            self.choose(kept)

    def choose(self, keys):
        self.session.select(keys)
        self.syncing = True
        self.table.selection_set(sorted(self.session.selected))
        if self.session.selected:
            self.table.see(sorted(self.session.selected)[0])
        self.syncing = False
        for key in self.session.keys:
            self.table.item(key, tags=("neighbor",) if key in self.session.neighbors else ())
        self.event_callback(
            dict(event="selection", residue_keys=sorted(self.session.selected), variant_id=self.session.variant_id)
        )
        self.draw()

    def table_selected(self, event):
        if not self.syncing:
            values = set(self.table.selection())
            if values != self.session.selected:
                self.choose(values)

    def variant_selected(self, event):
        self.session.select_variant(self.variant.get())
        self.choose(self.session.selected)
        row = self.session.document.predictions[
            self.session.document.predictions.variant_id.astype(str) == self.session.variant_id
        ]
        record = row.iloc[0]
        value = record.get("prediction", record.get("predicted", "unavailable"))
        try:
            number = float(value)
            value = f"{number:.4g}" if np.isfinite(number) else "Unavailable"
        except (TypeError, ValueError):
            value = "Unavailable"
        status = {"ok": "Complete", "base_only": "Base model only", "error": "Failed",
                  "pending": "Pending", "cancelled": "Cancelled"}.get(record.get("status"), "Unavailable")
        target = record.get("target", self.session.document.report.get("model", {}).get("target"))
        units = record.get("units") or {"pkd": "pKd", "ddg": "kcal/mol"}.get(target, "units unavailable")
        scope = self.session.document.report.get("model", {}).get("scope", record.get("model", ""))
        missing = [
            str(site.get("residue_key") or f"{site['chain']}:{site['position']}") + "=" + site["coverage"]
            for site in json.loads(record.sites)
            if site["coverage"] != "footprint"
        ]
        self.output.set(
            f"Variant {self.session.variant_id} · {record.get('mutations') or 'whole complex'} · "
            f"Prediction: {value} {units} · {status} · {scope}"
            + (" · Coverage: " + "; ".join(missing) if missing else "")
        )

    def search(self):
        q = self.query.get().strip().lower()
        self.choose([key for key in self.session.keys if q and q in key.lower()])

    def pick(self, event):
        keys = self.pick_maps.get(event.artist)
        if keys is not None and len(event.ind):
            key = str(keys[int(event.ind[0])])
            self.choose(self.session.selected ^ {key} if event.mouseevent.key == "shift" else [key])

    def pick_surface(self, event):
        if not event.dblclick or event.inaxes is not self.figure.axes[0] or event.xdata is None:
            return
        from topoppi.visualization.surface_rendering import ProjectedSurface

        triangles, labels = self.surface_pick
        if not len(triangles):
            return
        face, depth, weights = ProjectedSurface(triangles, event.inaxes.get_proj()).query(
            np.array([[event.xdata, event.ydata]])
        )
        if face[0] >= 0:
            key = str(labels[face[0], weights[0].argmax()])
            self.choose(self.session.selected ^ {key} if event.key == "shift" else [key])

    def colors(self, keys):
        return [
            "#CC6677"
            if k in self.session.selected
            else "#228833"
            if k in self.session.neighbors
            else "#4477AA"
            if self.session.side_by_key.get(k) == 0
            else "#AA8866"
            for k in keys
        ]

    def draw(self):
        from matplotlib.collections import PolyCollection
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        from topoppi.atlas.footprints import mesh_vertex_residue_labels
        from topoppi.visualization.footprint_rendering import footprint_geometry
        from topoppi.visualization.surface_rendering import surface_geometry

        self.session.condition = next(key for key, label in VIEW_LABELS.items() if label == self.condition.get())
        condition = self.session.condition
        camera = (
            (self.figure.axes[0].elev, self.figure.axes[0].azim)
            if self.figure.axes and hasattr(self.figure.axes[0], "elev")
            else None
        )
        limits = (
            [self.figure.axes[0].get_xlim(), self.figure.axes[0].get_ylim(), self.figure.axes[0].get_zlim()]
            if camera
            else None
        )
        self.figure.clear()
        self.pick_maps = {}
        ax = self.figure.add_subplot(121, projection="3d")
        triangles = []
        triangle_labels = []
        for partner in self.session.document.interface.partners:
            for patch in partner.patches:
                labels = mesh_vertex_residue_labels(patch, partner.residue_keys)
                geo = surface_geometry(patch, labels)
                triangles.append(geo["triangles"])
                triangle_labels.append(geo["corner_labels"])
                coll = Poly3DCollection(
                    geo["cells"], facecolors=self.colors(geo["labels"]), edgecolors="none", alpha=0.9
                )
                ax.add_collection3d(coll)
            # Residue centroids include sites outside the retained surface and
            # provide stable 3D picking without dependence on polygon z sorting.
            keys = sorted(set(map(str, partner.residue_keys)))
            xyz = np.array([partner.coords[np.asarray(partner.residue_keys) == k].mean(axis=0) for k in keys])
            dots = ax.scatter(*xyz.T, c=self.colors(keys), s=12, picker=5, depthshade=False)
            self.pick_maps[dots] = keys
        self.surface_pick = (
            (np.concatenate(triangles), np.concatenate(triangle_labels))
            if triangles
            else (np.empty((0, 3, 3)), np.empty((0, 3)))
        )
        coords = np.concatenate([p.coords for p in self.session.document.interface.partners])
        lo, hi = coords.min(axis=0), coords.max(axis=0)
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(hi - lo)
        ax.set_title("Native interface · select dots / double-click surface", fontsize=9)
        if camera:
            ax.view_init(*camera)
            ax.set_xlim(*limits[0])
            ax.set_ylim(*limits[1])
            ax.set_zlim(*limits[2])
        ax.set_xlabel("Å")
        ax.set_ylabel("Å")
        ax.set_zlabel("Å")
        right = self.figure.add_subplot(122)
        right.set_aspect("equal")
        right.axis("off")
        messages = []
        if condition == "G":
            keys = sorted(k for k in self.session.keys if self.session.contacts.get(k) or k in self.session.selected)
            # Deterministic bipartite barycentric ordering, same edges in every view.
            sides = [[k for k in keys if self.session.side_by_key[k] == s] for s in (0, 1)]
            for _ in range(8):
                for side in (0, 1):
                    rank = {k: i for i, k in enumerate(sides[1 - side])}
                    sides[side].sort(
                        key=lambda k: (
                            np.mean([rank[n] for n in self.session.contacts.get(k, {}) if n in rank])
                            if self.session.contacts.get(k)
                            else len(rank),
                            k,
                        )
                    )
            positions = {k: (s, i / max(1, len(side) - 1)) for s, side in enumerate(sides) for i, k in enumerate(side)}
            for a in sides[0]:
                for b in self.session.contacts.get(a, {}):
                    if b in positions:
                        strong = a in self.session.selected or b in self.session.selected
                        right.plot(
                            [0, 1],
                            [positions[a][1], positions[b][1]],
                            c="#228833" if strong else "#BBBBBB",
                            lw=1.5 if strong else 0.4,
                            alpha=1 if strong else 0.35,
                        )
            xy = np.array([positions[k] for k in keys])
            if len(xy):
                self.pick_maps[right.scatter(*xy.T, c=self.colors(keys), s=28, picker=5)] = keys
            for k in self.session.selected | self.session.neighbors:
                if k in positions:
                    right.text(*positions[k], k, fontsize=6)
            right.set_aspect("auto")
            right.set_title("Native residue-contact graph · ≤6 Å", fontsize=9)
        elif condition == "A":
            offset = 0.0
            for side in (0, 1):
                try:
                    key = (self.session.layout, side)
                    if key not in self.atlas_cache:
                        self.atlas_cache[key] = self.session.document.view(side=side, mode="atlas")
                    patches, viz, style = self.atlas_cache[key]
                except ValueError as exc:
                    messages.append(f"Partner {side}: {exc}")
                    continue
                for patch in patches:
                    labels = mesh_vertex_residue_labels(
                        patch, self.session.document.interface.partners[side].residue_keys
                    )
                    uv = patch.metadata.get("uv_global", patch.metadata.get("uv"))
                    geo = footprint_geometry(patch, uv, labels)
                    polygons = geo["polygons"].copy()
                    polygons[:, :, 0] += offset - polygons[:, :, 0].min()
                    offset = polygons[:, :, 0].max() + 4
                    artist = PolyCollection(
                        polygons, facecolors=self.colors(geo["labels"]), edgecolors="none", picker=True
                    )
                    right.add_collection(artist)
                    self.pick_maps[artist] = geo["labels"]
            right.autoscale_view()
            right.set_title("Surface atlas · all pieces linked by residue identity", fontsize=9)
        else:
            right.text(
                0.02,
                0.98,
                "Native contacts for selected residues\n\n"
                + "\n".join(
                    f"{k} → " + ", ".join(self.session.contacts.get(k, {})) for k in sorted(self.session.selected)
                ),
                va="top",
                fontsize=8,
                wrap=True,
                transform=right.transAxes,
            )
        self.status.set(
            f"{VIEW_LABELS[condition]} | Red: selected ({len(self.session.selected)}); green: contacting residues ({len(self.session.neighbors)}). Shift-click adds a residue. "
            + " ".join(messages)
        )
        self.context_text.configure(state="normal")
        self.context_text.delete("1.0", "end")
        lines = []
        table = self.session.document.residues
        for key in sorted(self.session.selected | self.session.neighbors):
            row = table[table.residue_key == key]
            sasa = row.iloc[0].get("buried_sasa_A2", "unavailable") if len(row) else "unavailable"
            lines.append(f"{key} · ΔSASA {sasa} Å² · contacts: " + ", ".join(self.session.contacts.get(key, {})))
        self.context_text.insert("1.0", "\n".join(lines))
        self.context_text.configure(state="disabled")
        self.figure.tight_layout()
        if self.expanded.get():
            right.set_visible(False)
            ax.set_position([0.04, 0.05, 0.92, 0.88])
        self.canvas.draw_idle()

    def show_measurements(self):
        window = tk.Toplevel(self.window)
        window.title("Experimental alanine-substitution measurements")
        text = tk.Text(window, width=110, height=25, wrap="word")
        text.pack(fill="both", expand=True)
        text.insert(
            "1.0",
            ("Experimental values apply to the listed measured substitutions.\n\n"
             + "\n\n".join(" | ".join(f"{k}: {v}" for k, v in r.items()) for r in self.measurements))
            if self.measurements else "Open an analysis with experimental measurements to view them here.",
        )
        text.configure(state="disabled")

    def save(self):
        path = filedialog.asksaveasfilename(parent=self.window, defaultextension=".json", initialfile="review.json")
        if path:
            self.session.save(path)
            self.event_callback(
                dict(
                    event="save_review",
                    path=path,
                    variant_id=self.session.variant_id,
                    residue_keys=sorted(self.session.selected),
                )
            )

    def open(self):
        path = filedialog.askopenfilename(parent=self.window, filetypes=[("Review", "*.json")])
        if path:
            try:
                reopened = ReviewSession.open(path)
                if reopened.document.report["analysis_id"] != self.session.document.report["analysis_id"]:
                    raise ValueError("Choose a review saved for the current analysis, or open its analysis.json in Affinity first.")
                if self.condition_locked:
                    reopened.condition = self.session.condition
                self.session = reopened
                self.condition.set(VIEW_LABELS[reopened.condition])
                self.variant.set(reopened.variant_id or "")
                self.choose(reopened.selected)
                if reopened.variant_id:
                    kept = set(reopened.selected)
                    self.variant_selected(None)
                    self.choose(kept)
                self.event_callback(
                    dict(
                        event="reopen_review",
                        path=path,
                        variant_id=reopened.variant_id,
                        residue_keys=sorted(reopened.selected),
                    )
                )
            except (ValueError, OSError) as exc:
                messagebox.showerror("Open review", str(exc), parent=self.window)

    def export_contacts(self):
        import pandas as pd

        path = filedialog.asksaveasfilename(parent=self.window, defaultextension=".csv", initialfile="contacts.csv")
        if path:
            pd.DataFrame(
                [
                    dict(source=a, target=b, atom_pair_count=n, cutoff_A=6.0)
                    for a in sorted(self.session.contacts)
                    if self.session.side_by_key[a] == 0
                    for b, n in self.session.contacts[a].items()
                ]
            ).to_csv(path, index=False)


def main(argv=None):
    import argparse

    from topoppi_affinity.workflow import load_analysis

    p = argparse.ArgumentParser()
    p.add_argument("analysis")
    p.add_argument("--view", choices=list(VIEW_CONDITIONS), default="atlas")
    a = p.parse_args(argv)
    path = a.analysis
    obj = json.loads(__import__("pathlib").Path(path).read_text()) if str(path).endswith(".json") else {}
    session = (
        ReviewSession.open(path)
        if obj.get("format") == "topoppi-affinity-review"
        else ReviewSession(load_analysis(path), condition=VIEW_CONDITIONS[a.view])
    )
    root = tk.Tk()
    root.withdraw()
    window = ReviewWindow(root, session)
    window.window.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()


if __name__ == "__main__":
    main()
