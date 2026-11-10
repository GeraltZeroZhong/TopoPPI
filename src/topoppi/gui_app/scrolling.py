"""Keep keyboard navigation visible inside the desktop settings panes."""

import tkinter as tk
from tkinter import ttk


def _inside(widget, content):
    while widget is not None:
        if widget is content:
            return True
        widget = getattr(widget, "master", None)
    return False


def reveal_focused_widget(canvas, content, widget):
    if not _inside(widget, content):
        return
    canvas.update_idletasks()
    bounds = canvas.bbox("all")
    if not bounds or bounds[3] <= bounds[1]:
        return
    top = widget.winfo_rooty() - canvas.winfo_rooty()
    bottom = top + widget.winfo_height()
    height = canvas.winfo_height()
    if top < 8 or widget.winfo_height() > height - 16:
        delta = top - 8
    elif bottom > height - 8:
        delta = bottom - height + 8
    else:
        return
    canvas.yview_moveto((canvas.canvasy(0) + delta - bounds[1]) / (bounds[3] - bounds[1]))


def scroll_settings_page(canvas, content, event, direction):
    if isinstance(event.widget, (tk.Text, tk.Listbox, ttk.Treeview, ttk.Combobox)):
        return None
    if _inside(event.widget, content):
        canvas.yview_scroll(direction, "pages")
        return "break"
    return None
