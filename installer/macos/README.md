# macOS application

TopoPPI 2.1 provides disk images for Apple Silicon and Intel Macs:

```text
TopoPPI-2.1-macos-arm64.dmg       Apple Silicon
TopoPPI-2.1-macos-x86_64.dmg      Intel
```

Download an image from the [v2.1 release](https://github.com/GeraltZeroZhong/TopoPPI/releases/tag/v2.1),
or follow the [local build instructions](#build-locally).

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
~/Library/Application Support/TopoPPI/2.1-<architecture>
```

Later launches reuse that runtime and start more quickly. The disk image already
contains Python, the scientific dependencies and OptCuts. Runtime preparation
and analysis use these local files and work offline.

TopoPPI, OptCuts, and third-party license notices are stored inside the app at
`Contents/Resources`. The machine-readable inventory records packages whose
upstream metadata leaves the license file unspecified, for review during
release preparation.

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
4. Move the `2.1-arm64` or `2.1-x86_64` folder to the Trash.
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

Use the package version directly from the source tree:

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

## Release contract

The `macOS App` GitHub Actions workflow builds on native Apple Silicon and Intel
macOS 15 runners. Each job runs a full command-line mapping with its native
OptCuts executable, creates the application disk image, expands the packaged
runtime, and creates a footprint map with the packaged executable. It then
reopens the saved atlas through the CLI and exercises Tk highlighting, numerical
annotation, figure export and atlas saving. The image passes to the central
`Publish` workflow after these checks.

The packaged launcher and OptCuts executable target macOS 12 or later. The app
uses ad-hoc code signing, with the first-open steps described above.
