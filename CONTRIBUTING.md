# Contributing to TopoPPI

This guide covers TopoPPI 3.0 development from a source checkout to a tested
change. The [v2.1 tag](https://github.com/GeraltZeroZhong/TopoPPI/tree/v2.1) identifies
the TopoPPI paper version and retains its code and documentation.

## Set up the development environment

TopoPPI uses Python 3.10. Create the Conda environment and install the project in editable mode:

```bash
conda env create -f environment.yml
conda activate topoppi-dev
python -m pip install -e ".[dev,benchmark,meshio,representation,analysis]"
```

On Linux x86-64, install the bundled OptCuts executable before running the complete mapping pipeline:

```bash
bash tools/OptCuts/install_optcuts.sh
command -v OptCuts_bin
```

Windows and macOS development requires a native OptCuts build. The instructions live beside their installers:

- [Windows installer guide](./installer/windows/README.md)
- [macOS application guide](./installer/macos/README.md)

## Work on the affinity plugin

Install the optional plugin from the repository root. The full test suite also
requires PyTorch and PyTorch Geometric for its neural-model tests:

```bash
python -m pip install -e . -e "./plugins/affinity[dev]"
python -m pip install "torch>=2.2" "torch-geometric>=2.5"
python -m pytest plugins/affinity/tests
```

The bundled affinity and mutation-effect models run on CPU with the plugin's
standard dependencies. Install `./plugins/affinity[models]` when working on
optional neural or sequence models. The [plugin test guide](plugins/affinity/tests/README.md) describes synthetic fixtures and optional dependencies. The [workflow guide](plugins/affinity/docs/workflow.md) covers prediction from structures and saved-result review.

## Work on the hotspot plugin

```bash
python -m pip install -e ./plugins/affinity -e ./plugins/hotspot
python -m pytest plugins/hotspot/tests
```

## Work on representations and analysis

`src/topoppi_representation/` contains representation-specific data objects,
native-mesh, UV and direct-table codecs and physical descriptors. It uses Core
atoms, partner identities and meshes. Install
`.[representation]` for these interfaces. The [representation guide](docs/representation.md)
and [analysis script guide](scripts/analysis/README.md) document optional readout
dependencies and explicit input/partition formats.

Keep scientific feature definitions, schema ordering and encoded-byte accounting
stable when reorganizing code. Name task-specific chemistry and normalization
schemes explicitly. Packages use installed resources, and scripts accept
user-selected inputs and output directories.

## Test the standalone platform boundary

Use a fresh Python 3.10 environment with `pip install .`, then run
`python -I tests/core_install_smoke.py`. Install `./examples/plugin` and rerun
with `--plugin` to verify discovery and execution without editing Core. This
smoke test checks the default Chemistry plugin, ProLIF detection and the absence
of optional prediction packages. CI separately provisions geometry dependencies,
installs the wheel with `--no-deps`, and runs with `--without-chemistry-deps` to
check geometry-only operation. Use the regular development environment above
for the complete regression suite.

Editable installation registers the plugin and backend entry-point metadata
needed for discovery; `PYTHONPATH` controls module imports. Chemistry code lives
in `src/topoppi_chemistry`. Its [backend contract](docs/chemistry.md#adding-a-backend)
supports future perception implementations without modifying Core.

New extensions should follow the [plugin contract](docs/architecture.md#plugin-contract).
Keep registration modules lightweight and pass shared Core objects to capabilities.

## Make a focused change

- Keep code, comments, commit messages, and repository documentation in English.
- Preserve scientific behavior unless the change intentionally updates it and includes matching tests and documentation.
- Prefer direct code paths and clear errors. Add validation where it protects a real user-facing boundary.
- Keep software fixtures small, synthetic, and redistributable. Retain documented teaching structures, product screenshots and the three default model packages. Store research datasets, intermediate weights and run outputs separately from the source tree.
- Update the changelog and user documentation when a command, default, output, or installation step changes.

## Run the checks

Run the same core checks used in continuous integration:

```bash
python -m pytest
python -m ruff check src tests installer
```

On a headless Linux machine, run the tests with `xvfb-run -a python -m pytest`
after installing Xvfb and xauth. CI uses this command so the Tk editing tests
execute alongside the pipeline tests. Windows and macOS run Tk tests through
their native display systems.

Tests that launch the external OptCuts executable use the `requires_optcuts` marker. Dataset-scale or long-running tests use `slow`. A new test should carry one of these markers only when its runtime or dependency warrants it.

For packaging changes, also build and inspect both distributions:

```bash
python -m build
python -m twine check dist/*
```

For user-facing command changes, exercise the relevant help page and one representative workflow. Examples include:

```bash
topoppi --help
topoppi render --help
topoppi-benchmark --help
topoppi-install-optcuts --help
```

For the three default task models, install the Core and both plugin wheels into
a clean environment. Copy `installer/smoke_affinity.py` and
`installer/smoke_binding.py` together outside the checkout, then run:

```bash
python smoke_binding.py --no-gui
```

This executes actual affinity, mutation-effect and hotspot inference on an
independently constructed peptide pair, checks finite outputs and units, and
reopens saved results after the source file is removed. Python network requests
are blocked during the smoke. Use `--require-gui` with a display to also require
Tk window creation. `smoke_affinity.py --neural` separately checks an explicitly
installed optional neural stack without loading weights.

The platform workflows also run `installer/smoke_footprints.py` against a saved
atlas from the installed CLI. This checks the packaged GUI's atlas opening,
highlighting, CSV annotation and save operations.

## Open a pull request

Describe the user problem, the chosen behavior, and the checks you ran. Include screenshots for visible desktop changes and a small output example for schema or command-line changes. Keep unrelated cleanup in a separate change so reviewers can trace each result to its implementation.

## GUI embedding

The application class is `topoppi.gui_app.InterfaceApp`. Integrations importing the former `ProtSurfApp` class should update that import and constructor name. CLI entry points remain unchanged.
