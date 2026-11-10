# Affinity tests

Install `./plugins/affinity[dev]` for the test tools, POT surface-transport solvers and CatBoost regressors. Install Torch and PyTorch Geometric to include the neural-model tests. From the repository root, run `python -m pytest plugins/affinity/tests`. Public unit tests use deterministic numerical arrays, analytically constructed meshes and the independently constructed `tiny_complex.pdb` coordinates. Random inputs use explicit seeds.

Tests also use the documented 1BVK structure to check complete multi-chain selection and example export. Numerical tests construct their fitted models from synthetic inputs; default-model tests use the bundled weights. Tests with explicit optional-dependency checks skip when their software is unavailable; GUI tests require a display, for example `xvfb-run -a python -m pytest plugins/affinity/tests`.
