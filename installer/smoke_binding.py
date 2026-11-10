"""Exercise all three installed default models on an independent synthetic PDB.

Copy this file and smoke_affinity.py together outside the checkout, then run it
with the Python interpreter that has the Core and both plugin wheels installed.
The default Hotspot workflow requires an available native OptCuts executable.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from smoke_affinity import check_affinity, offline


def check_hotspot(directory):
    import numpy as np
    from topoppi_hotspot.gui_workflow import RankingSettings, load_ranking, run_ranking

    from topoppi.examples import copy_example

    directory = Path(directory)
    source = copy_example(directory / "input", synthetic=True)
    settings = RankingSettings(str(source), "A", "B", "A", "", str(directory / "ranking"))
    document = run_ranking(settings)
    assert len(document.rows) == 1
    assert document.metadata["model_id"] == "topoppi-hotspot-mutation-surface"
    assert document.metadata["target"] == "ddg"
    assert document.metadata["units"] == "kcal/mol"
    assert document.metadata["temperature_K"] == 298.15
    assert np.isfinite([float(row["rank_score"]) for row in document.rows]).all()
    np.testing.assert_allclose([float(row["rank_score"]) for row in document.rows],
                               [float(row["predicted_ddg_kcal_mol"]) for row in document.rows])
    source.unlink()
    restored = load_ranking(settings.output_dir)
    assert restored.geometry is not None
    assert restored.rows == [{key: str(value) for key, value in row.items()} for row in document.rows]
    patches, visualizer, style = restored.view(restored.rows[0]["residue_key"])
    assert patches and visualizer is not None and style["view"] == "surface"
    restored.export(directory / "ranking-values.csv", annotations=True)
    assert (directory / "ranking-values.csv").read_text().startswith("residue,value\n")
    assert "torch" not in sys.modules, "Default binding analysis must run without PyTorch"
    return {"rows": len(document.rows), "saved_geometry": True, "annotation_export": True}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-gui", action="store_true", help="Require Tk window creation")
    parser.add_argument("--no-gui", action="store_true", help="Run all scientific workflows without a Tk window")
    args = parser.parse_args(argv)
    with tempfile.TemporaryDirectory(prefix="topoppi-binding-smoke-") as folder, offline():
        result = check_affinity(Path(folder) / "affinity", require_gui=args.require_gui, gui=not args.no_gui)
        result["tasks"]["hotspot"] = check_hotspot(Path(folder) / "hotspot")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
