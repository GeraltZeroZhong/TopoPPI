import tkinter as tk
from copy import deepcopy
from tkinter import colorchooser, messagebox

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.collections import PathCollection, PolyCollection

from topoppi.errors import ConfigurationError
from topoppi.visualization.visualizer import select_patches_for_display


class PlotMixin:
    def redraw_plot(self):
        if self._busy:
            self.log("The current run will use its saved display settings.")
            return
        if self._pending_single_run is not None:
            self._render_pending_result()
            return
        successful_run = self._successful_single_run
        if successful_run is None:
            return
        self._capture_surface_camera()
        try:
            style = self.get_style_config()
        except (ValueError, ConfigurationError) as exc:
            self.log(f"Could not update style: {exc}")
            return
        self.log("Updating plot style...")
        self.update_plot(
            successful_run["viz"],
            successful_run.get("all_patches", successful_run["patches"]),
            style,
            complete_task=False,
        )

    def _render_pending_result(self, style=None):
        pending = self._pending_single_run
        try:
            style = dict(style) if style is not None else self.get_style_config()
            style["min_points"] = pending["params"].get("min_points", pending["viz"].config.min_points)
            for key in ("label_offsets", "marker_color_overrides", "residue_color_overrides"):
                style[key] = dict(pending["style"].get(key, {}))
        except (ValueError, ConfigurationError) as exc:
            self._finish_task()
            self.log(f"Calculation complete. Adjust the display settings and apply the style: {exc}")
            return False
        return self.update_plot(
            pending["viz"], pending["patches"], style, complete_task=True,
            run_params=pending["params"], run_manifest=pending["manifest"],
            all_patches=pending["all_patches"],
        )

    def on_pick(self, event):
        if self._surface_view() is not None:
            return
        artist = event.artist
        if not isinstance(artist, (PathCollection, PolyCollection)) or self._busy or self._pending_single_run is not None:
            return
        gid = artist.get_gid()
        successful_run = self._successful_single_run
        viz = successful_run["viz"] if successful_run else None
        if gid and viz and gid in viz.artist_map:
            if isinstance(artist, PolyCollection) and viz.last_style.get("annotation_values") is not None:
                self.log("Numeric annotations control region colors. Clear annotations to recolor residues.")
                return
            color = colorchooser.askcolor(title=f"Color for {gid}")[1]
            if color:
                self._mark_style_custom()
                target_objs = viz.artist_map[gid]
                if isinstance(artist, PolyCollection):
                    residue_key = target_objs["residue_key"]
                    self.residue_color_overrides[residue_key] = color
                    for objects in viz.artist_map.values():
                        if objects.get("residue_key") == residue_key and objects.get("collection") is not None:
                            objects["collection"].set_facecolor(color)
                else:
                    self.marker_color_overrides[gid] = color
                    target_objs["scatter"].set_facecolor(color)
                self.current_canvas.draw()
                self._remember_interactive_style()

    def _remember_interactive_style(self):
        self._capture_surface_camera()
        successful_run = self._successful_single_run
        style = dict(successful_run["viz"].last_style)
        style.update(
            {
                "label_offsets": dict(self.label_offsets),
                "marker_color_overrides": dict(self.marker_color_overrides),
                "residue_color_overrides": dict(self.residue_color_overrides),
            }
        )
        successful_run["style"] = style
        successful_run["viz"].last_style = dict(style)

    @staticmethod
    def _label_anchor(objects):
        anchor = objects.get("anchor")
        if anchor is not None:
            return anchor
        scatter = objects.get("scatter")
        return scatter.get_offsets()[0] if scatter is not None else None

    def on_mouse_press(self, event):
        successful_run = self._successful_single_run
        if self._busy or self._pending_single_run is not None or successful_run is None or not event.inaxes:
            return
        surface = self._surface_view()
        if surface is not None:
            self._drag_state = None
            if getattr(event, "dblclick", False) and event.button == 1:
                self._recolor_surface_residue(surface, event)
            return
        for gid, objs in successful_run["viz"].artist_map.items():
            txt = objs.get("text")
            if txt is None:
                continue
            contains, _ = txt.contains(event)
            if contains:
                self._drag_state = {"gid": gid}
                break

    def on_mouse_move(self, event):
        if self._surface_view() is not None:
            return
        if not self._drag_state or not event.inaxes or event.xdata is None or event.ydata is None:
            return
        gid = self._drag_state["gid"]
        objs = self._successful_single_run["viz"].artist_map.get(gid, {})
        txt = objs.get("text")
        if txt is None:
            return
        txt.set_position((event.xdata, event.ydata))
        connector = objs.get("connector")
        pt = self._label_anchor(objs)
        if connector is not None and pt is not None:
            connector.set_data([pt[0], event.xdata], [pt[1], event.ydata])
        self.current_canvas.draw_idle()

    def on_mouse_release(self, event):
        if self._surface_view() is not None:
            self._drag_state = None
            self._capture_surface_camera()
            return
        if not self._drag_state:
            return
        gid = self._drag_state["gid"]
        objs = self._successful_single_run["viz"].artist_map.get(gid, {})
        txt = objs.get("text")
        pt = self._label_anchor(objs)
        if txt is not None and pt is not None:
            tx, ty = txt.get_position()
            self.label_offsets[gid] = (float(tx - pt[0]), float(ty - pt[1]))
            self._remember_interactive_style()
        self._drag_state = None

    def _surface_view(self):
        return getattr(getattr(self, "current_fig", None), "_topoppi_surface", None)

    def _capture_surface_camera(self):
        surface = self._surface_view()
        successful_run = self._successful_single_run
        if surface is None or successful_run is None or self._pending_single_run is not None:
            return None
        camera = deepcopy(surface["capture_camera"]())
        successful_run["style"]["surface_camera"] = camera
        successful_run["viz"].last_style["surface_camera"] = deepcopy(camera)
        self.loaded_atlas_style["surface_camera"] = deepcopy(camera)
        return camera

    def reset_surface_view(self):
        if self._busy or self._pending_single_run is not None:
            return
        surface = self._surface_view()
        if surface is None:
            self.loaded_atlas_style.pop("surface_camera", None)
            return
        surface["reset_camera"]()
        self._capture_surface_camera()
        self.current_canvas.draw_idle()

    def _recolor_surface_residue(self, surface, event):
        viz = self._successful_single_run["viz"]
        if viz.last_style.get("annotation_values") is not None:
            self.log("Numeric annotations control region colors. Clear annotations to recolor residues.")
            return
        gid = surface["pick_residue"](event)
        if gid is None or gid not in viz.artist_map:
            return
        residue = viz.artist_map[gid]["residue_key"]
        color = colorchooser.askcolor(title=f"Color for {residue}")[1]
        if color:
            self._mark_style_custom()
            self.residue_color_overrides[residue] = color
            residue_gids = [key for key, objects in viz.artist_map.items() if objects["residue_key"] == residue]
            self.marker_color_overrides = {
                key: value for key, value in self.marker_color_overrides.items()
                if not any(key == uid or key.startswith(uid + "__piece_") for uid in residue_gids)
            }
            self._remember_interactive_style()
            self.redraw_plot()

    def update_plot(
        self,
        viz,
        patches,
        style,
        complete_task=False,
        *,
        run_params=None,
        run_manifest=None,
        all_patches=None,
    ):
        previous_render = {name: getattr(viz, name, None) for name in
                           ("artist_map", "last_style", "last_report", "capture_surface_camera")}
        try:
            all_patches = list(patches if all_patches is None else all_patches)
            params = run_params if complete_task else self._successful_single_run["params"]
            style = dict(style)
            style["min_points"] = params.get("min_points", viz.config.min_points)
            patches, _counts = select_patches_for_display(
                all_patches, viz, map_style=style.get("map_style", "markers"),
                min_points=params.get("min_points", viz.config.min_points),
                view=style.get("view", "atlas"),
            )
            fig = viz.plot_patches(patches, show=False, style_config=style)
        except Exception as exc:
            for name, value in previous_render.items():
                setattr(viz, name, value)
            if complete_task:
                pending = getattr(self, "_pending_single_run", None)
                message = (
                    f"Calculation complete. The atlas is ready to save.\n"
                    f"Adjust the display settings and click Apply Style to draw it.\n\n{exc}"
                    if pending is not None else f"Failed to generate plot: {exc}"
                )
                self.show_error(self._previous_result_message(message))
            else:
                self.log(f"Failed to update plot style: {exc}")
            return False

        annotation = dict(viz.last_report)

        old_fig = self.current_fig
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()
        if old_fig is not None and old_fig is not fig:
            plt.close(old_fig)

        self.current_fig = fig
        self.current_canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        self.current_canvas.draw()
        toolbar_frame = tk.Frame(self.canvas_frame, bg="#ffffff")
        toolbar_frame.pack(side=tk.TOP, fill=tk.X, padx=8, pady=(4, 0))
        hidden_toolbar = tk.Frame(self.canvas_frame)
        self.current_toolbar = NavigationToolbar2Tk(self.current_canvas, hidden_toolbar, pack_toolbar=False)
        self.current_toolbar.update()
        self._build_plot_toolbar(toolbar_frame)
        self.current_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.current_canvas.mpl_connect("pick_event", self.on_pick)
        self.current_canvas.mpl_connect("button_press_event", self.on_mouse_press)
        self.current_canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.current_canvas.mpl_connect("button_release_event", self.on_mouse_release)
        if complete_task:
            manifest = dict(run_manifest)
            manifest["visualization_annotation"] = annotation
            successful_run = {
                "figure": fig,
                "manifest": manifest,
                "params": dict(run_params),
                "patches": patches,
                "all_patches": all_patches,
                "log": list(self.current_run_log),
                "style": dict(viz.last_style),
                "viz": viz,
            }
            self._successful_single_run = successful_run
            self._pending_single_run = None
            self._run_style = None
            self._restore_atlas_style(viz.last_style)
            self.var_annotation_target.set("Current map")
            if self._next_run_annotations is not None and self._next_run_annotations["context"] == self._result_context(successful_run):
                self._next_run_annotations = None
            self._finish_task()
            if successful_run["params"].get("auto_save"):
                try:
                    self._auto_save_current_figure(successful_run)
                except Exception as exc:
                    self.log(f"Auto-save failed: {exc}")
        else:
            successful_run = self._successful_single_run
            successful_run["figure"] = fig
            successful_run["patches"] = patches
            successful_run["all_patches"] = all_patches
            successful_run["manifest"]["visualization_annotation"] = annotation
            successful_run["style"] = dict(viz.last_style)
            self._refresh_action_states()
        annotation = viz.last_report
        interaction_source = str(annotation.get("interaction_residue_source", "interaction"))
        source_label = {"prolif": "ProLIF", "geometric": "geometric", "geometric_fallback": "geometric"}.get(
            interaction_source, "recorded"
        )
        residue_count = int(annotation.get("displayed_residue_count", 0))
        self.log(
            f"Showing {len(patches)} {'patch' if len(patches) == 1 else 'patches'} and "
            f"{residue_count} {'residue' if residue_count == 1 else 'residues'} "
            f"({int(annotation.get('patch_interaction_residue_count', 0))} with {source_label} contacts)."
        )
        return True

    def _build_plot_toolbar(self, parent):
        surface = self._surface_view() is not None
        actions = [
            ("Reset view", self.reset_surface_view) if surface else ("Home", self.current_toolbar.home),
            ("Back", self.current_toolbar.back),
            ("Forward", self.current_toolbar.forward),
            ("Pan", self.current_toolbar.pan),
            ("Zoom", self.current_toolbar.zoom),
        ]
        controls = tk.Frame(parent, bg="#ffffff")
        controls.pack(anchor=tk.W)
        for label, command in actions:
            tk.Button(controls, text=label, command=command, relief=tk.FLAT, padx=8, pady=2).pack(
                side=tk.LEFT, padx=(0, 4)
            )
        hint = ("Drag to rotate; double-click a residue to change its color." if surface
                else "Drag labels to move them; click residues to change colors.")
        tk.Label(parent, text=hint,
                 bg="#ffffff", fg="#52616B", font=(self.config.font_family, 9)).pack(anchor=tk.W)

    def _close_current_figure(self):
        if self.current_fig is not None:
            plt.close(self.current_fig)
            self.current_fig = None
        self.current_canvas = None
        self.current_toolbar = None

    def show_error(self, msg):
        self._restore_display_highlights()
        self._finish_task()
        self.log(f"Error: {msg}")
        messagebox.showerror("Pipeline Error", msg)
