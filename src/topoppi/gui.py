"""GUI entry point for TopoPPI."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from topoppi import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="topoppi-gui",
        description="Create, annotate, and export protein interface maps in the TopoPPI desktop app.",
        epilog="Choose a PDB/mmCIF structure to create a map, or use Open Atlas to continue editing a saved map.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--structure", help="Preselect a PDB or mmCIF structure.")
    parser.add_argument("-A", "--chain-a", help="Preselect surface/first partner chain(s).")
    parser.add_argument("-B", "--chain-b", help="Preselect second partner chain(s).")
    parser.add_argument("--output-dir", help="Preselect the output directory.")
    parser.add_argument("--optcuts-bin", help="Use this OptCuts executable.")
    parser.add_argument("--plugin", help="Open an installed plugin interface after startup.")
    parser.add_argument("--plugin-setup", help=argparse.SUPPRESS)
    parser.add_argument("--affinity", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--affinity-setup", help=argparse.SUPPRESS)
    return parser


def _set_window_icon(root) -> None:
    """Set the Tk window icon when packaged assets are available."""
    import sys
    import tkinter as tk
    from importlib import resources

    try:
        asset_files = resources.files("topoppi.assets")
        if sys.platform.startswith("win"):
            with resources.as_file(asset_files / "topoppi.ico") as icon_path:
                root.iconbitmap(default=str(icon_path))

        with resources.as_file(asset_files / "topoppi.png") as icon_path:
            icon_image = tk.PhotoImage(file=str(icon_path))
        root.iconphoto(True, icon_image)
        root._topoppi_icon_image = icon_image
    except (FileNotFoundError, ModuleNotFoundError, OSError, tk.TclError):
        return


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    import tkinter as tk

    import matplotlib

    from topoppi.config import DEFAULT_GUI_CONFIG

    matplotlib.use("TkAgg", force=True)

    from topoppi.gui_app import InterfaceApp

    root = tk.Tk(className="TopoPPI")
    _set_window_icon(root)
    root.tk.call("tk", "scaling", DEFAULT_GUI_CONFIG.tk_scaling)
    app = InterfaceApp(root, config=DEFAULT_GUI_CONFIG)
    if args.structure:
        app.var_input_path.set(args.structure)
        app._populate_chain_choices(args.structure)
    for value, variable in ((args.chain_a, app.var_chain_a), (args.chain_b, app.var_chain_b),
                            (args.output_dir, app.var_output_dir)):
        if value is not None:
            variable.set(value)
    if args.optcuts_bin:
        app.entry_optcuts_bin.delete(0, tk.END)
        app.entry_optcuts_bin.insert(0, args.optcuts_bin)
    app._form_changed()
    from topoppi.gui_compat import extension_arguments

    extension, setup_file = extension_arguments(args)
    if extension:
        def open_extension():
            import json
            from pathlib import Path
            from tkinter import messagebox

            path = Path(setup_file) if setup_file else None
            try:
                setup = json.loads(path.read_text(encoding="utf-8")) if path else None
                app.open_extension(extension, setup=setup)
                if path:
                    path.unlink()
            except (OSError, ValueError, RuntimeError) as exc:
                retained = f"\nThe saved input is retained at {path}" if path else ""
                messagebox.showerror("Open plugin interface", f"{exc}{retained}", parent=root)
        root.after_idle(open_extension)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
