"""Local ranking and saved results used by the desktop interface."""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class RankingSettings:
    structure: str
    chains_a: str
    chains_b: str
    partner: str
    model_dir: str
    output_dir: str
    embedding_registry: str = ""


@dataclass
class RankingDocument:
    directory: Path
    metadata: dict
    rows: list[dict]
    geometry: object = None

    def export(self, path, *, annotations=False):
        """Export the displayed results or a Core-compatible residue/value CSV."""
        fields = ["residue", "value"] if annotations else list(self.rows[0]) if self.rows else ["residue_key", "rank_score"]
        with Path(path).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(({"residue": row["residue_key"], "value": row["rank_score"]} for row in self.rows)
                             if annotations else self.rows)

    def view(self, selected=None):
        if self.geometry is None:
            raise ValueError("Open a ranking folder containing interface.npz to display its saved surface.")
        if not self.rows:
            raise ValueError("Choose a partner with standard non-alanine interface residues to display ranking colors.")
        group = tuple(self.metadata["selected_partner"])
        side = next(index for index, partner in enumerate(self.geometry.partners) if partner.chains == group)
        patches = self.geometry.partners[side].patches
        if not patches:
            raise ValueError("This partner has an empty surface interface. The residue table retains its ranking results.")
        selected_keys = {selected} if isinstance(selected, str) else set(selected or ())
        highlights = tuple(row["residue_key"] for row in self.rows
                           if row["residue_key"] in selected_keys and float(row.get("footprint_area_A2", 0)) > 0)
        return patches, self.geometry.visualizer(side), {
            "map_style": "footprints", "view": "surface", "residue_scope": "patch",
            "footprint_labels": "highlighted", "highlight_residues": highlights,
            "annotation_values": {row["residue_key"]: float(row["rank_score"]) for row in self.rows},
            "annotation_label": ("Predicted Ala ΔΔG (kcal/mol)" if self.metadata.get("target") == "ddg"
                                 else "Hotspot ranking score"), "show_seams": False,
        }


def load_ranking(path):
    """Read saved ranking outputs; model files and source coordinates are optional."""
    from topoppi import load_interface

    directory = Path(path).expanduser().resolve()
    if directory.is_file():
        directory = directory.parent
    metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    with (directory / "rankings.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not {"residue_key", "rank_score", "rank_min", "rank_max"}.issubset(reader.fieldnames or []):
            raise ValueError("Choose the rankings.csv and metadata.json folder produced by Hotspot.")
        rows = list(reader)
    expected_rows = metadata.get("candidate_coverage", {}).get("n_ranked_candidates")
    if not rows or (expected_rows is not None and len(rows) != expected_rows):
        raise ValueError("Saved rankings table is incomplete. Restore rankings.csv from the complete results folder.")
    native = directory / metadata.get("native_interface", "interface.npz")
    geometry = load_interface(native) if native.is_file() else None
    return RankingDocument(directory, metadata, rows, geometry)


def run_ranking(settings, *, prepared=None, context=None):
    """Run the installed capability on compatible Core geometry, using CPU."""
    import topoppi
    from topoppi.interface import chain_group

    from .inference import model_requirements

    context = context or topoppi.ExecutionContext()
    context.check_cancelled()
    from .resources import resolve_model_dir

    model = resolve_model_dir(settings.model_dir)
    source = Path(settings.structure).expanduser().resolve()
    groups = (chain_group(settings.chains_a), chain_group(settings.chains_b))
    selected = chain_group(settings.partner)
    if set(groups[0]) & set(groups[1]):
        raise ValueError("Assign each chain to one partner group.")
    if selected not in groups:
        raise ValueError("Choose one complete partner group to rank.")
    output = Path(settings.output_dir).expanduser().resolve()
    if any((output / filename).exists() for filename in (
        "rankings.csv", "metadata.json", "interface.npz", "candidate_coverage.csv", "features.csv", "representation.bin",
    )):
        raise ValueError("Choose a new results folder to preserve the saved ranking.")
    requirements = model_requirements(model)
    if requirements["needs_embeddings"] and not settings.embedding_registry:
        raise ValueError("This model uses sequence embeddings. Select its prepared embedding registry in Sequence model options.")
    if settings.embedding_registry and not Path(settings.embedding_registry).expanduser().is_dir():
        raise ValueError("Choose an existing embedding registry folder.")
    try:
        output.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ValueError(f"Choose a writable results folder. The selected destination could not be prepared: {exc}") from exc
    protocol = requirements["feature_protocol"]
    same_input = (prepared is not None and Path(prepared.source["path"]).resolve() == source
                  and {partner.chains for partner in prepared.partners} == set(groups))
    compatible = (same_input and asdict(prepared.surface_config) == protocol["surface"]
                  and asdict(prepared.topology_config) == protocol["topology"]
                  and all(partner.mesh is not None for partner in prepared.partners))
    if compatible:
        geometry = prepared
        context.emit("geometry", message="Reusing the current interface geometry.")
    else:
        context.emit("geometry", message="Preparing both partner surfaces with the model settings.")
        if same_input:
            interface = prepared.interface
        else:
            if not settings.structure or not source.is_file():
                raise ValueError("Choose a PDB or mmCIF structure to rank.")
            interface = topoppi.get_interface(topoppi.load_structure(source), *groups)
        geometry = topoppi.compute_geometry(interface, surface=topoppi.SurfaceConfig(**protocol["surface"]),
            topology=topoppi.TopologyConfig(**protocol["topology"]), cancel_event=context.cancel_event)
    context.check_cancelled()
    context.emit("ranking", message="Calculating residue features and ranking the selected partner.")
    frame, metadata = topoppi.run("hotspot.rank", geometry, partner=selected, model_dir=model, output_dir=output,
        embedding_registry=settings.embedding_registry or None, device="cpu", context=context)
    context.check_cancelled()
    topoppi.save_interface(geometry, output / "interface.npz")
    metadata["native_interface"] = "interface.npz"
    metadata["geometry_reused"] = bool(compatible)
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return RankingDocument(output, metadata, frame.to_dict(orient="records"), geometry)
