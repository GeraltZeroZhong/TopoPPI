"""Exercise installed Core/Affinity defaults outside the source checkout."""
from __future__ import annotations

import argparse
import json
import socket
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace


@contextmanager
def offline():
    """Make accidental Python network requests fail during the local smoke."""
    original = socket.socket.connect

    def blocked(*args, **kwargs):
        raise AssertionError("Installed default prediction attempted network access")

    socket.socket.connect = blocked
    try:
        yield
    finally:
        socket.socket.connect = original


def check_affinity(directory, *, require_gui=False, gui=True):
    import numpy as np
    import topoppi_affinity
    from topoppi_affinity.workflow import analyze, load_analysis

    import topoppi
    from topoppi.examples import copy_example

    assert topoppi.__version__ == topoppi_affinity.__version__
    assert "torch" not in sys.modules, "Default imports must leave neural runtimes unloaded"
    directory = Path(directory)
    structure = copy_example(directory / "input", synthetic=True)
    geometry = topoppi.compute_geometry(topoppi.get_interface(topoppi.load_structure(structure), "A", "B"))
    outputs = {}
    for task, mutations, column in (("absolute", (), "predicted_pkd"),
                                     ("mutation", ("A:F2A",), "predicted_ddg_kcal_mol")):
        output = directory / task
        report = analyze(structure, "A", "B", task=task, mutations=mutations,
                         output=output, prepared=geometry, log=lambda *_: None)
        assert report["prediction_status"] == "ok", report.get("prediction_error", report)
        document = load_analysis(output / "analysis.json")
        assert document.interface is not None
        assert len(document.predictions) == 1
        assert document.predictions.status.eq("ok").all()
        assert np.isfinite(document.predictions[column].to_numpy(float)).all()
        patches, visualizer, style = document.view(side=0, mode="surface")
        assert patches and style["view"] == "surface" and visualizer is not None
        assert (output / "predictions.csv").is_file()
        outputs[task] = {"rows": len(document.predictions), "saved_geometry": True, "unit_column": column}
    structure.unlink()
    for task in outputs:
        assert load_analysis(directory / task / "analysis.json").interface is not None
    assert "torch" not in sys.modules, "Bundled physical models must run without PyTorch"
    gui_status = "not_requested"
    if gui:
        import tkinter as tk

        from topoppi_affinity.gui import open_panel

        try:
            root = tk.Tk()
        except tk.TclError:
            if require_gui:
                raise
            gui_status = "unavailable"
        else:
            root.withdraw()
            panel = open_panel(SimpleNamespace(root=root))
            try:
                root.update_idletasks()
                assert "torch" not in sys.modules
                gui_status = "passed"
            finally:
                panel.close()
                root.destroy()
    return {"version": topoppi.__version__, "tasks": outputs, "gui": gui_status}


def check_neural_imports():
    """Optional advanced-model environment check, without loading weights."""
    import numpy as np
    import torch
    from esm.models.esmc import ESMC  # noqa: F401

    np.testing.assert_array_equal(torch.from_numpy(np.array([1., 2.], dtype=np.float32)).numpy(), [1., 2.])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-gui", action="store_true", help="Fail when a Tk display is unavailable")
    parser.add_argument("--no-gui", action="store_true", help="Run prediction and saved-result checks only")
    parser.add_argument("--neural", action="store_true", help="Also check an explicitly installed optional neural stack")
    args = parser.parse_args(argv)
    with tempfile.TemporaryDirectory(prefix="topoppi-affinity-smoke-") as folder, offline():
        result = check_affinity(folder, require_gui=args.require_gui, gui=not args.no_gui)
        if args.neural:
            check_neural_imports()
            result["optional_neural_imports"] = "passed"
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
