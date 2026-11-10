"""Portable selection and native-contact correspondence, independent of layout."""

from __future__ import annotations

import json
import os
from pathlib import Path

from topoppi.atlas.footprints import source_atom_residue_labels

from topoppi_affinity.workflow import load_analysis

VIEW_CONDITIONS = {"structure": "D", "contacts": "G", "atlas": "A"}
VIEW_LABELS = {"D": "Structure", "G": "Contact graph", "A": "Interface atlas"}


class ReviewSession:
    def __init__(self, document, *, condition="A", layout="primary"):
        if condition not in ("D", "G", "A"):
            raise ValueError("Choose the structure, contact graph or atlas view")
        if document.interface is None:
            raise ValueError("Review requires a saved native interface")
        if layout not in ("primary", "residue_aware"):
            raise ValueError("Unknown saved atlas configuration")
        if layout == "residue_aware":
            records = json.loads((document.path.parent / "secondary_layout.json").read_text())
            document.report["layout"] = {
                str(r["side"]): {
                    "status": r["status"],
                    "atlas": r.get("atlas"),
                    "error": r.get("error", "Saved atlas unavailable"),
                }
                for r in records
            }
        self.layout = layout
        self.document = document
        self.condition = condition
        self.selected = set()
        self.variant_id = None
        self.annotation = "footprint_area_A2"
        self.side = 0
        self.keys = set()
        self.side_by_key = {}
        self.contacts = {}
        for side, partner in enumerate(document.interface.partners):
            for key in source_atom_residue_labels(partner.atoms):
                self.keys.add(str(key))
                self.side_by_key[str(key)] = side
            self.contacts.update(document.interface.contacts(side))

    def select(self, keys):
        keys = set(map(str, keys))
        if keys - self.keys:
            raise ValueError(f"Unknown residues: {sorted(keys - self.keys)}")
        self.selected = keys
        return self.neighbors

    def select_variant(self, variant_id):
        rows = self.document.predictions
        rows = rows[rows.variant_id.astype(str) == str(variant_id)]
        if len(rows) != 1:
            raise ValueError(f"Unknown or ambiguous variant ID: {variant_id}")
        sites = json.loads(rows.iloc[0].sites)
        self.select(s["residue_key"] for s in sites if s.get("residue_key") in self.keys)
        self.variant_id = str(variant_id)

    @property
    def neighbors(self):
        return {k for selected in self.selected for k in self.contacts.get(selected, {})}

    def save(self, path):
        path = Path(path).resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        report = self.document.report
        obj = dict(
            format="topoppi-affinity-review",
            version=1,
            analysis=os.path.relpath(self.document.path, path.parent),
            analysis_id=report["analysis_id"],
            selected_residue_keys=sorted(self.selected),
            selected_variant_id=self.variant_id,
            view=next(view for view, condition in VIEW_CONDITIONS.items() if condition == self.condition),
            layout=self.layout,
            side=self.side,
            annotation=self.annotation,
            contact_cutoff_A=6.0,
            units=report.get("units"),
            model=report.get("model"),
        )
        path.write_text(json.dumps(obj, indent=2) + "\n")
        return path

    @classmethod
    def open(cls, path):
        path = Path(path).resolve()
        obj = json.loads(path.read_text())
        if obj.get("format") != "topoppi-affinity-review" or obj.get("version") != 1:
            raise ValueError("Unsupported review document")
        session = cls(
            load_analysis(path.parent / obj["analysis"]),
            condition=VIEW_CONDITIONS[obj["view"]] if "view" in obj else obj["condition"],
            layout=obj.get("layout", "primary"),
        )
        if session.document.report["analysis_id"] != obj["analysis_id"]:
            raise ValueError("Review references a different analysis")
        if obj["selected_variant_id"] is not None:
            session.select_variant(obj["selected_variant_id"])
        session.select(obj["selected_residue_keys"])
        session.side = obj["side"]
        session.annotation = obj["annotation"]
        return session
