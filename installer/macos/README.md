# macOS application

This guide describes TopoPPI 3.0.0 source builds and installation. Build the current
checkout using the [local instructions](#build-locally). Disk-image filenames
use the source version and target architecture:

```text
TopoPPI-<version>-macos-arm64.dmg       Apple Silicon
TopoPPI-<version>-macos-x86_64.dmg      Intel
```

The TopoPPI paper uses **v2.1**. Its [release](https://github.com/GeraltZeroZhong/TopoPPI/releases/tag/v2.1)
and [version-specific macOS instructions](https://github.com/GeraltZeroZhong/TopoPPI/blob/v2.1/installer/macos/README.md)
remain available. Select the matching release when reproducing that version.

Open **Apple menu > About This Mac** when you need to check the processor. A Mac
with an Apple chip uses `arm64`; a Mac showing an Intel processor uses `x86_64`.

## Install and start

1. Open the disk image that matches the Mac.
2. Drag **TopoPPI** into **Applications**.
3. Try to open **TopoPPI**. If macOS blocks it, open **System Settings >
   Privacy & Security**, choose **Open Anyway** for TopoPPI, and confirm
   **Open**. Older macOS releases may also offer **Open** through the app's
   Control-click menu. These builds use ad-hoc signing.
4. Keep the preparation window open while TopoPPI expands its bundled runtime.
   This commonly takes several minutes on the first launch.

The app stores its prepared runtime at:

```text
~/Library/Application Support/TopoPPI/<version>-<architecture>
```

Later launches reuse that runtime and start more quickly. The disk image already
contains Python, the scientific dependencies and OptCuts. Runtime preparation
and analysis use these local files and work offline.

TopoPPI, OptCuts, and third-party license notices are stored inside the app at
`Contents/Resources`. A machine-readable inventory provides the bundled package and license metadata.

## Binding analysis in 3.0

Choose **Binding analysis...** for **Binding affinity**, **Mutation effects**, or
**Rank hotspots**. Affinity supplies the first two tasks; Hotspot supplies the
third. Install matching plugin wheels through **Plugins → Manage Plugins... →
Install from file...**, or include the plugins in the runtime before building
the app. Each task has an official bundled CPU model for offline prediction
from local structures.

The [desktop guide](../../docs/desktop.md#binding-analysis) covers the three
workflows and residue-to-variant selection. The [Affinity guide](../../plugins/affinity/docs/workflow.md)
and [Hotspot guide](../../plugins/hotspot/README.md) describe outputs, model scope
and optional custom models. Advanced neural models use their separately
installed dependencies and declared resources.

## Upgrade

Quit TopoPPI, open the disk image for the new release, and replace the existing
app in **Applications**. The first launch of a new version prepares a separate
versioned runtime. Older runtime folders remain available until you remove them.

Figures, manifests, structures, and benchmark results saved in other folders
remain in place.

## Repair startup

Startup details are written to:

```text
~/Library/Logs/TopoPPI/launcher.log
```

When the app reports a startup failure:

1. Quit TopoPPI.
2. In Finder, choose **Go > Go to Folder**.
3. Open `~/Library/Application Support/TopoPPI`.
4. Move the folder matching the affected app version and architecture to the Trash.
5. Open TopoPPI again and keep the preparation window open.

This rebuilds the packaged runtime from the copy inside the application.

## Uninstall

Quit TopoPPI and move `/Applications/TopoPPI.app` to the Trash. To remove the
prepared runtimes as well, open `~/Library/Application Support/TopoPPI` in
Finder and move that folder to the Trash.

Analysis output remains in the folders you selected. The recent-file list is
stored at `~/.topoppi/gui_recent.json`; remove the `.topoppi` folder when you
also want to clear that preference data. The launcher log can be removed from
`~/Library/Logs/TopoPPI`.

## Build locally

Build on the target architecture in a Conda or Micromamba environment containing
TopoPPI and its runtime dependencies. Install these build tools as well:

- Xcode Command Line Tools
- Git
- CMake
- `conda-pack`
- network access for the pinned OptCuts source checkout

To include the three prediction tasks, install `./plugins/affinity` and
`./plugins/hotspot` in that environment before packaging. Use the package version
directly from the source tree:

```bash
python -m pip install "prolif>=2.0"
python -m pip install --no-deps .
version="$(PYTHONPATH=src python -c 'from topoppi._version import __version__; print(__version__)')"
architecture="$(uname -m)"

bash tools/OptCuts/build_residue_aware_optcuts.sh \
  "release-assets/OptCuts_bin-macos-${architecture}"
bash installer/macos/build_app.sh \
  "$version" \
  "$architecture" \
  "$CONDA_PREFIX" \
  "release-assets/OptCuts_bin-macos-${architecture}" \
  release-assets
```

The disk image is written under `release-assets/`. Open it and complete one
first-launch check on the same architecture before sharing it.

When testing another build with the same application version, remove its
prepared runtime as described in [Repair startup](#repair-startup). The next
launch expands the runtime from the updated application.

## Packaged components

The disk image contains the application, a Python runtime, scientific dependencies and the matching native OptCuts executable. Apple Silicon and Intel images are built separately. The launcher and executable target macOS 12 or later, and the app uses ad-hoc signing with the first-open steps above.

The `macOS App` workflow builds both architectures and exercises mapping, saved-atlas editing and figure export from the packaged runtime before passing the images to the `Publish` workflow.
