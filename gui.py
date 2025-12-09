import tkinter as tk
from tkinter import filedialog, messagebox, ttk, colorchooser
import threading
import sys
import os
import shutil
import subprocess
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Ensure src is in python path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.io_loader import PDBLoader
from src.surface import SurfaceGenerator
from src.topology import TopologyManager
from src.parameterization import Parameterizer
from src.visualizer import InterfaceVisualizer

class ProtSurfApp:
    def __init__(self, root):
        self.root = root
        self.root.title("ProtSurf2D - Protein Interface Mapper")
        self.root.geometry("1200x950") 

        self.cached_viz = None
        self.cached_patches = None
        self._picking = False 

        # Interaction Types for Arpeggio
        self.interaction_types_list = [
            'salt_bridge',
            'cation_pi',
            'pi_stack',
            'hbond',
            'sulfur_pi',
            'hydrophobic',
            'halogen_bond',
            'metal_complex'
        ]
        
        self.default_active = {
            'salt_bridge', 'hbond', 'hydrophobic', 'cation_pi', 'pi_stack'
        }
        self.interaction_vars = {}

        self.paned_window = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.left_frame = ttk.Frame(self.paned_window, width=340)
        self.right_frame = ttk.Frame(self.paned_window)
        self.paned_window.add(self.left_frame, weight=1)
        self.paned_window.add(self.right_frame, weight=4)

        self._init_controls()
        self._init_plot_area()
        self._init_status_bar()

    def _init_controls(self):
        # 1. Input
        lbl_frame = ttk.LabelFrame(self.left_frame, text="1. Input Data", padding=10)
        lbl_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # PDB Input
        ttk.Label(lbl_frame, text="PDB Structure:").pack(anchor=tk.W)
        self.entry_file = ttk.Entry(lbl_frame)
        self.entry_file.pack(fill=tk.X, pady=2)
        ttk.Button(lbl_frame, text="Browse PDB...", command=self.browse_file).pack(anchor=tk.E)

        # Arpeggio Input
        ttk.Label(lbl_frame, text="Arpeggio JSON (Auto-generated if empty):").pack(anchor=tk.W, pady=(5,0))
        self.entry_arpeggio = ttk.Entry(lbl_frame)
        self.entry_arpeggio.pack(fill=tk.X, pady=2)
        ttk.Button(lbl_frame, text="Browse JSON...", command=self.browse_arpeggio).pack(anchor=tk.E)

        # 2. Parameters
        param_frame = ttk.LabelFrame(self.left_frame, text="2. Analysis Parameters", padding=10)
        param_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(param_frame, text="Chain A (Surf):").grid(row=0, column=0, sticky='w')
        self.entry_chain_a = ttk.Entry(param_frame, width=10)
        self.entry_chain_a.insert(0, "A")
        self.entry_chain_a.grid(row=0, column=1, pady=2)

        ttk.Label(param_frame, text="Chain B (Lig):").grid(row=1, column=0, sticky='w')
        self.entry_chain_b = ttk.Entry(param_frame, width=10)
        self.entry_chain_b.insert(0, "B")
        self.entry_chain_b.grid(row=1, column=1, pady=2)

        ttk.Label(param_frame, text="Cutoff (Å):").grid(row=2, column=0, sticky='w')
        self.entry_cutoff = ttk.Entry(param_frame, width=10)
        self.entry_cutoff.insert(0, "5.0")
        self.entry_cutoff.grid(row=2, column=1, pady=2)

        ttk.Label(param_frame, text="Grid Res (Å):").grid(row=3, column=0, sticky='w')
        self.entry_res = ttk.Entry(param_frame, width=10)
        self.entry_res.insert(0, "1.0")
        self.entry_res.grid(row=3, column=1, pady=2)

        ttk.Label(param_frame, text="Sigma:").grid(row=4, column=0, sticky='w')
        self.entry_sigma = ttk.Entry(param_frame, width=10)
        self.entry_sigma.insert(0, "1.5")
        self.entry_sigma.grid(row=4, column=1, pady=2)

        # 3. Style Controls
        style_frame = ttk.LabelFrame(self.left_frame, text="3. Visualization Style", padding=10)
        style_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.var_color_type = tk.BooleanVar(value=True) 
        self.chk_type = ttk.Checkbutton(style_frame, text="Color by Interaction Type", 
                                        variable=self.var_color_type, command=self.toggle_color_mode)
        self.chk_type.pack(anchor=tk.W, pady=2)

        # Interaction Type Filter
        self.filter_frame = ttk.LabelFrame(style_frame, text="Show Interactions", padding=5)
        self.filter_frame.pack(fill=tk.X, pady=5)
        
        for i, itype in enumerate(self.interaction_types_list):
            var = tk.BooleanVar(value=(itype in self.default_active))
            self.interaction_vars[itype] = var
            chk = ttk.Checkbutton(self.filter_frame, text=itype, variable=var, command=self.redraw_plot)
            chk.grid(row=i//2, column=i%2, sticky='w', padx=2)

        f_frame = ttk.Frame(style_frame)
        f_frame.pack(fill=tk.X, pady=5)
        ttk.Label(f_frame, text="Font:").pack(side=tk.LEFT)
        self.combo_font = ttk.Combobox(f_frame, values=["Arial", "Times New Roman", "Courier New", "sans-serif"], width=10)
        self.combo_font.current(3)
        self.combo_font.pack(side=tk.LEFT, padx=2)
        self.spin_size = ttk.Spinbox(f_frame, from_=5, to=20, width=4)
        self.spin_size.set(9)
        self.spin_size.pack(side=tk.LEFT)

        self.residue_color = "#ff0000"
        self.btn_color = tk.Button(f_frame, text="Def. Color", bg=self.residue_color, command=self.choose_color, relief=tk.RAISED, state=tk.DISABLED, width=10)
        self.btn_color.pack(side=tk.RIGHT, padx=5)

        btn_frame = ttk.Frame(self.left_frame)
        btn_frame.pack(fill=tk.X, padx=10, pady=10)
        self.btn_run = ttk.Button(btn_frame, text="Run Full Analysis", command=self.start_analysis)
        self.btn_run.pack(fill=tk.X, pady=5)
        self.btn_redraw = ttk.Button(btn_frame, text="Update Style Only", command=self.redraw_plot, state=tk.DISABLED)
        self.btn_redraw.pack(fill=tk.X, pady=5)
        self.progress = ttk.Progressbar(self.left_frame, mode='indeterminate')
        self.progress.pack(fill=tk.X, padx=10, pady=5)

    def _init_plot_area(self):
        self.canvas_frame = ttk.Frame(self.right_frame)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)
        self.current_canvas = None
        lbl = ttk.Label(self.canvas_frame, text="Load a PDB then click Run.\nArpeggio will be invoked automatically if JSON is missing.", font=("Arial", 12))
        lbl.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

    def _init_status_bar(self):
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        self.status = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

    def log(self, message):
        self.status_var.set(message)
        print(f"[GUI Log] {message}")
        self.root.update_idletasks()

    def browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("PDB Files", "*.pdb"), ("CIF Files", "*.cif"), ("All Files", "*.*")])
        if filename:
            self.entry_file.delete(0, tk.END)
            self.entry_file.insert(0, filename)
            # Auto-clear Arpeggio field if user picks new PDB, to force re-calculation or manual re-select
            self.entry_arpeggio.delete(0, tk.END)

    def browse_arpeggio(self):
        filename = filedialog.askopenfilename(filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")])
        if filename:
            self.entry_arpeggio.delete(0, tk.END)
            self.entry_arpeggio.insert(0, filename)

    def choose_color(self):
        color = colorchooser.askcolor(color=self.residue_color, title="Select Residue Color")[1]
        if color:
            self.residue_color = color
            self.btn_color.config(bg=color)
            
    def toggle_color_mode(self):
        if self.var_color_type.get():
            self.btn_color.config(state=tk.DISABLED)
            for child in self.filter_frame.winfo_children():
                child.configure(state='normal')
        else:
            self.btn_color.config(state=tk.NORMAL)
            for child in self.filter_frame.winfo_children():
                child.configure(state='disabled')
        if self.cached_viz:
            self.redraw_plot()

    def get_style_config(self):
        active_types = [t for t, var in self.interaction_vars.items() if var.get()]
        return {
            'color': self.residue_color,
            'font_family': self.combo_font.get(),
            'font_size': int(self.spin_size.get()),
            'color_by_type': self.var_color_type.get(),
            'active_types': active_types
        }

    def start_analysis(self):
        pdb_path = self.entry_file.get()
        if not pdb_path or not os.path.exists(pdb_path):
            messagebox.showerror("Error", "Please select a valid PDB file.")
            return
            
        params = {
            'path': pdb_path,
            'chain_a': self.entry_chain_a.get().strip(),
            'chain_b': self.entry_chain_b.get().strip(),
            'arpeggio': self.entry_arpeggio.get().strip(),
            'cutoff': float(self.entry_cutoff.get()),
            'res': float(self.entry_res.get()),
            'sigma': float(self.entry_sigma.get())
        }
        self.btn_run.config(state=tk.DISABLED)
        self.btn_redraw.config(state=tk.DISABLED)
        self.progress.start(10)
        self.log("Starting analysis pipeline...")
        threading.Thread(target=self.run_pipeline, args=(params,), daemon=True).start()

    def redraw_plot(self):
        if not self.cached_viz or not self.cached_patches: return
        style = self.get_style_config()
        self.log("Updating plot style...")
        self.update_plot(self.cached_viz, self.cached_patches, style)

    def generate_arpeggio_interactions(self, pdb_path):
        """
        Attempts to run pdbe-arpeggio to generate interaction data.
        Runs directly in the current environment.
        """
        self.log("Arpeggio JSON missing. Attempting to generate...")
        
        pdb_dir = os.path.dirname(os.path.abspath(pdb_path))
        pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
        
        # Predicted output file (pdbe-arpeggio typically outputs {name}.json)
        expected_json = os.path.join(pdb_dir, f"{pdb_name}.json")
        
        # If it already exists from a previous run, use it
        if os.path.exists(expected_json):
            self.log(f"Found existing Arpeggio output: {os.path.basename(expected_json)}")
            return expected_json

        # Command construction
        # Assume pdbe-arpeggio is installed in the current PATH (via pip in environment.yml)
        cmd = ["pdbe-arpeggio", pdb_path]
        
        success = False
        
        try:
            self.log("Running: pdbe-arpeggio ...")
            # Run without 'conda run', using current env context
            subprocess.run(cmd, check=True, cwd=pdb_dir, capture_output=True)
            success = True
        except subprocess.CalledProcessError as e:
            self.log(f"Arpeggio generation failed: {e}")
            print(f"Subprocess Error (stderr): {e.stderr.decode()}")
        except FileNotFoundError:
            self.log("Error: 'pdbe-arpeggio' command not found. Did you install dependencies?")

        if success and os.path.exists(expected_json):
            self.log("Arpeggio JSON generated successfully.")
            
            # Update the GUI entry field for clarity (thread-safe way)
            self.root.after(0, lambda: self.entry_arpeggio.insert(0, expected_json))
            return expected_json
        else:
            self.log("Could not generate Arpeggio JSON. Falling back to geometric heuristics.")
            return None

    def run_pipeline(self, params):
        try:
            # --- Step 0: Arpeggio Handling ---
            arpeggio_file = params.get('arpeggio')
            
            # If path is empty or invalid, try to generate it
            if not arpeggio_file or not os.path.exists(arpeggio_file):
                generated_json = self.generate_arpeggio_interactions(params['path'])
                if generated_json:
                    arpeggio_file = generated_json
                else:
                    arpeggio_file = None # Explicitly None to trigger heuristic fallback
            
            # --- Step 1: Loading ---
            self.log("Loading PDB structure...")
            loader = PDBLoader(params['path'])
            coords_A, atoms_A = loader.get_chain_data(params['chain_a'])
            coords_B, atoms_B = loader.get_chain_data(params['chain_b'])

            # --- Step 2: Surface ---
            self.log("Generating molecular surface...")
            surf_gen = SurfaceGenerator(coords_A)
            mesh_A = surf_gen.generate_mesh(grid_resolution=params['res'], sigma=params['sigma'])
            if mesh_A is None: raise ValueError("Surface generation failed.")

            # --- Step 3: Topology ---
            self.log("Extracting interface patches...")
            topo = TopologyManager(mesh_A, coords_B)
            patches = topo.get_interface_patches(distance_cutoff=params['cutoff'])
            if not patches: raise ValueError("No interface found with current cutoff.")

            # --- Step 4: Parameterization ---
            self.log(f"Flattening {len(patches)} patches...")
            param = Parameterizer()
            valid_patches = []
            for p in patches:
                uv = param.flatten_patch(p)
                if uv is not None:
                    p.metadata['uv'] = uv
                    valid_patches.append(p)
            
            if not valid_patches: raise ValueError("LSCM Parameterization failed for all patches.")

            # --- Step 5: Visualization ---
            self.log("Rendering visualization...")
            
            viz = InterfaceVisualizer(
                chain_A_atoms=atoms_A, 
                chain_A_coords=coords_A, 
                chain_B_coords=coords_B, 
                chain_B_atoms=atoms_B,
                chain_a_id=params['chain_a'],
                chain_b_id=params['chain_b'],
                arpeggio_file=arpeggio_file
            )
            
            self.cached_viz = viz
            self.cached_patches = valid_patches
            self.root.after(0, lambda: self.finish_success())
        except Exception as e:
            err_msg = str(e)
            self.root.after(0, lambda: self.show_error(err_msg))

    def finish_success(self):
        style = self.get_style_config()
        self.update_plot(self.cached_viz, self.cached_patches, style)
        self.btn_redraw.config(state=tk.NORMAL)

    def on_pick(self, event):
        if self._picking: return
        artist = event.artist
        gid = artist.get_gid()
        if gid and self.cached_viz and gid in self.cached_viz.artist_map:
            self._picking = True
            try:
                color = colorchooser.askcolor(title=f"Color for {gid}")[1]
                if color:
                    target_objs = self.cached_viz.artist_map[gid]
                    target_objs['scatter'].set_facecolor(color)
                    target_objs['scatter'].set_edgecolor(color)
                    self.current_canvas.draw()
            finally:
                self._picking = False

    def update_plot(self, viz, patches, style):
        for widget in self.canvas_frame.winfo_children(): widget.destroy()
        fig = viz.plot_patches(patches, show=False, style_config=style)
        if fig is None:
            self.show_error("Failed to generate plot.")
            return
        self.current_canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        self.current_canvas.draw()
        self.current_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.current_canvas.mpl_connect('pick_event', self.on_pick)
        self.progress.stop()
        self.btn_run.config(state=tk.NORMAL)
        self.log(f"Success! Displaying {len(patches)} patches.")

    def show_error(self, msg):
        self.progress.stop()
        self.btn_run.config(state=tk.NORMAL)
        self.log("Error occurred.")
        messagebox.showerror("Pipeline Error", msg)

if __name__ == "__main__":
    root = tk.Tk()
    app = ProtSurfApp(root)
    root.mainloop()
