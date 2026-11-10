"""Open the map's chemical annotation controls through the plugin interface."""


class ChemistryController:
    def __init__(self, host):
        self.host = host

    def open(self, setup=None):
        from topoppi.gui_app.scrolling import reveal_focused_widget

        self.host.var_settings_page.set("advanced")
        self.host._sync_settings_page()
        control = self.host.combo_interaction_source
        control.focus_set()
        self.host.root.update_idletasks()
        reveal_focused_widget(self.host.sidebar_canvas, self.host.left_frame, control)
        self.host.status_var.set("Choose an interaction backend for your next interface map.")
        return control


def create_controller(host):
    return ChemistryController(host)
