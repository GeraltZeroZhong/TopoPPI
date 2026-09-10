"""Exercise 2D and 3D editing with the Python runtime installed by a release build."""

from __future__ import annotations

import argparse
import tkinter as tk
from pathlib import Path
from unittest.mock import patch

import numpy as np

from topoppi.atlas.uv import as_corner_uv
from topoppi.gui_app import ProtSurfApp
from topoppi.visualization.atlas_io import load_atlas


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("atlas", type=Path, help="Atlas produced by the installed CLI")
    parser.add_argument("output_dir", type=Path, help="Directory for the edited figure and atlas")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = load_atlas(args.atlas)
    image_path = args.output_dir / "desktop-footprints.svg"
    atlas_path = args.output_dir / "desktop-footprints.npz"
    surface_image = args.output_dir / "desktop-surface.png"
    surface_atlas = args.output_dir / "desktop-surface.npz"
    values_path = args.output_dir / "values.csv"

    root = tk.Tk()
    root.withdraw()
    app = ProtSurfApp(root)
    try:
        with (
            patch("topoppi.gui_app.ui_mixin.messagebox.showinfo"),
            patch("topoppi.gui_app.ui_mixin.messagebox.showerror", side_effect=AssertionError),
        ):
            with patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=str(args.atlas)):
                app.open_atlas()
            root.update()
            assert app._successful_single_run is not None, "The installed GUI did not open the atlas."
            visualizer = app._successful_single_run["viz"]
            assert visualizer.last_style["map_style"] == "footprints"
            key = next(iter(visualizer.artist_map.values()))["residue_key"]
            app.var_highlight_residues.set(key)
            app.redraw_plot()
            values_path.write_text(f"residue,value\n{key},1.2\n", encoding="utf-8")
            with patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=str(values_path)):
                app.browse_annotations()
            root.update()
            assert len(app.current_fig.axes) == 2, "The numeric view has no colorbar."
            with patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=str(image_path)):
                app.save_figure()
            with patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=str(atlas_path)):
                app.save_atlas()
            app.combo_view.set("3D interface")
            app._view_changed()
            root.update()
            axis = app.current_fig._topoppi_surface["axis"]
            axis.view_init(elev=48, azim=-62)
            app.current_canvas.draw()
            with patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=str(surface_image)):
                app.save_figure()
            with patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=str(surface_atlas)):
                app.save_atlas()

        edited = load_atlas(atlas_path)
        assert image_path.stat().st_size > 0
        assert edited.style["highlight_residues"] == [key]
        assert edited.style["annotation_values"] == {key: 1.2}
        assert len(source.patches) == len(edited.patches)
        for original, reopened in zip(source.patches, edited.patches, strict=True):
            np.testing.assert_array_equal(as_corner_uv(original), as_corner_uv(reopened))
        surface = load_atlas(surface_atlas)
        assert surface_image.stat().st_size > 0
        assert surface.style["view"] == "surface"
        assert surface.style["surface_camera"]["elevation"] == 48
        assert surface.style["surface_camera"]["azimuth"] == -62
        assert surface.style["annotation_values"] == {key: 1.2}
        for original, reopened in zip(source.patches, surface.patches, strict=True):
            np.testing.assert_array_equal(original.vertices, reopened.vertices)
            np.testing.assert_array_equal(as_corner_uv(original), as_corner_uv(reopened))
        print(f"Installed desktop footprint editing passed: {image_path}")
        print(f"Installed desktop 3D editing passed: {surface_image}")
    finally:
        app.close()


if __name__ == "__main__":
    main()
