# Changelog

All notable changes to TopoPPI are documented here.

## [2.1] - 2026-09-10

- Added a native 3D interface view in the CLI, GUI and Python API, with illuminated surface geometry and optional triangular mesh, residue borders and optimized seams.
- Shared residue markers, footprints, highlights and numerical annotations between the 2D atlas and 3D surface views. Every retained patch is shown at its original relative position in 3D.
- Added GUI camera rotation, pan and zoom, orthographic and perspective projection, a reset control, and double-click recoloring of visible residues.
- Saved the selected view, camera and mesh settings in editable atlas files, with computation-free switching and export through the GUI and `topoppi render`.
- Added CLI view, projection, camera and mesh options, plus explicit controls to restore hidden seams and residue borders when restyling saved atlases.
- Added a README gallery of four native output styles, current desktop screenshots and direct installation and usage links.

## [2.0] - 2026-09-09

- Added native Residue footprints mode in CLI and GUI, with exact per-corner residue regions, complete patch coverage, selected-residue highlighting, labels, boundaries and seams.
- Added CSV numerical annotations, missing-value coloring, configurable value scales, and shared SVG/PDF/PNG/TIFF exports.
- Added self-contained atlas NPZ export, GUI Save Atlas / Open Atlas, and `topoppi render` for styling saved coordinates without rerunning the solver.
- Added GUI region recoloring and label dragging, with numerical annotation colors kept consistent with their value scale.
- Added explicit geometric interaction selection for runs using distance-defined optimization weights.
- Updated the Windows launchers to resolve their environment from the installation directory, including paths with spaces and Unicode characters. Same-version setup reinstalls the application package for repair.
- Added a persistent Windows installation log and propagated environment setup errors to the installer's exit status.
- Supported Unicode Windows installation directories with a temporary cache drive, automatic mapping cleanup, and UTF-8 command output.
- Enabled the existing Tk interaction tests in Linux CI and added installed-runtime footprint editing and atlas rendering to platform smoke runs.
- Kept ProLIF calculations for each chain pair in the calling process, avoiding nested process pools in GUI and benchmark workers.
- Preserved completed computation when a display setting needs correction, allowing style recovery and atlas export without repeating optimization.
- Preserved residue highlights, manual colors and label positions after failed or cancelled runs, and waited for active tasks to finish cleanup when closing the desktop app.
- Kept residue annotation scope consistent across CLI rendering and GUI atlas editing.
- Separated the displayed map's input and annotations from the next run, captured settings when a run starts, and adapted desktop controls to smaller screens.
- Added UTF-8 BOM support for annotation CSV files and colorbar extensions for values outside selected limits, with corresponding residue counts in the rendering report.
- Restored marker interaction thresholds when rendering saved atlases while retaining every optimized patch for footprint mode.
- Expanded CLI help for saved-atlas rendering and footprint controls, and synchronized the public guides with the 2.0 release.

## [1.3] - 2026-08-30

### Desktop and installation

- Added self-contained macOS disk images for Apple Silicon and Intel, including the matching native OptCuts executable and an offline runtime.
- Added visible first-launch preparation and recovery guidance to the macOS app.
- Improved the Windows setup flow with a console-free GUI launcher, a persistent configured command prompt, an automation wrapper, direct installation of the bundled OptCuts executable, complete ProLIF dependency checks, and reliable bootstrap failure reporting.
- Added native platform smoke runs to the Windows and macOS release workflows.
- Expanded platform documentation for installation, upgrades, repair, removal, and local release builds.
- Added the installed version, project help, issue reporting, and an output-folder control to the Basic desktop workflow.
- Bound saving and redraw operations to the last successful single-run snapshot, so cancelled, failed, and benchmark tasks cannot mix figures with another run's settings.

### Interface mapping

- Added the residue-footprint-aware OptCuts objective as the complete TopoPPI method, identified in benchmark output as `residue_aware_optcuts`, together with a matched geometry-only ablation.
- Added seam-preserving per-face-corner UV coordinates across OBJ input and output, metrics, atlas packing, and visualization.
- Added LSCM initialization for OptCuts, automatic-initialization comparison, timeout and cancellation support, binary and artifact hashes, and audited upstream provenance.
- Added deterministic area-matched chart packing with recorded transforms, spacing, overlap, utilization, and waste measurements.
- Updated distortion metrics to use a common similarity scale, original 3-D area weighting, Jacobian singular values, and global-reflection-corrected flip counts.
- Added internal cut-edge seam counts and normalized 3-D seam length.
- Made supplied and generated ProLIF records the shared source for interaction membership, residue weights, GUI patch filtering, labels, and colors. Geometric interaction inference remains available through an explicit option.
- Bound formal benchmark residue weights to the manifest's checked ProLIF record across comparative and operational profiles.
- Added chemically perceived, explicitly hydrogenated RDKit copies for ProLIF fingerprinting while preserving the source structure coordinates.
- Added automatic ProLIF generation for PDB, CIF, and mmCIF inputs, including multi-character mmCIF chain identifiers.
- Aligned CLI, GUI, Python, and Windows application output around the same interaction and visualization rules.
- Unified the GUI, CLI, and Python interface cutoff at `4.0 Å` and changed the default display threshold to one interaction residue for small interfaces.
- Moved structure, chain, and output validation ahead of ProLIF and OptCuts work; missing-chain errors now list available protein chains.
- Placed generated ProLIF records beside the selected output image and added `--version` to each public command.

### Benchmarking and reproducibility

- Rebuilt comparisons around one provenance-tracked 3-D face domain shared by LSCM, harmonic, SLIM, spherical, cylindrical, and the selected OptCuts methods.
- Added isolated repetitions, process-tree memory sampling, fixed thread settings, warm-ups, repeatability checks, and explicit resource metadata.
- Added formal manifests, coordinate audits, input and interaction checksums, deterministic chain selection, inclusion logs, per-patch retention, pLDDT propagation, per-face samples, and compressed provenance mappings.
- Added cluster-aware paired inference, cluster bootstrap confidence intervals, effect sizes, Benjamini-Hochberg correction, and JSON `null` handling for unavailable values.
- Added surface-grid and memory preflight, controlled adaptive grid spacing, read-only benchmark preflight, fingerprint-bound resume, and one-factor or factorial sensitivity plans.
- Added experimental and AlphaFold DB paired-structure preparation, confidence and geometry strata, and paired transport analysis tools.
- Added a manifest preparation command that generates ProLIF records in bulk and fills the formal `prolif_file` and `prolif_sha256` fields.
- Added concise terminal summaries, field-specific configuration errors, `topoppi-benchmark verify`, safe output-root checks for every profile, and a runnable quick-start configuration.

### Maintenance

- Removed obsolete compatibility paths, repeated validation, unreachable guards, unused wrappers, and duplicate file helpers across the core pipeline, GUI, and publication tools.
- Expanded regression and real-binary integration coverage for mapping, seams, packing, provenance, benchmark evidence, installers, and end-to-end execution.
- Reworked the README, benchmark schema, platform guides, OptCuts notices, contribution guide, and publication-tool help around task-oriented user workflows.

## [1.2] - 2026-05-09

- Added Windows x86-64 bootstrap installer scaffolding for one-click TopoPPI installation.
- Added Windows-aware OptCuts artifact resolution in `topoppi-install-optcuts`.
- Added Linux OptCuts runtime sidecar installation for `libigl_stb_image.so`.
- Added release workflow support for optional Windows OptCuts artifacts and Windows setup executables.
- Documented Windows installer requirements, release artifacts, and troubleshooting paths.

## [1.1] - 2026-05-09

- Improved GUI mode switching, sticky run controls, styling, run logs, validation, figure auto-save, manifests, and Matplotlib navigation.
- Added benchmark `resume`/`new`/`overwrite` modes, worker selection, config fingerprints, and invalid-input reports.
- Made benchmark preprocessing honor configured chain IDs consistently.
- Added ProLIF metadata validation and richer single-run reproducibility manifests.
- Hardened release packaging, PyPI metadata, and OptCuts installation checks.
- Added the `topoppi-install-optcuts` downloader for Linux x86-64 GitHub release artifacts.

## [1.0] - 2026-04-06

- Initial public release of TopoPPI.
