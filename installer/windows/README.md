# Windows installer

This guide describes TopoPPI 3.0.0 source builds and installation. Build the current
checkout using the [local instructions](#build-locally); the setup filename uses
the version stored in the source tree:

```text
TopoPPI-<version>-windows-x86_64-setup.exe
```

The TopoPPI paper uses **v2.1**. Its [release](https://github.com/GeraltZeroZhong/TopoPPI/releases/tag/v2.1)
and [version-specific Windows instructions](https://github.com/GeraltZeroZhong/TopoPPI/blob/v2.1/installer/windows/README.md)
remain available. Select the matching release when reproducing that version.

The setup program creates a dedicated Python 3.10 environment under
`%LOCALAPPDATA%\TopoPPI`, installs TopoPPI and its scientific dependencies,
copies the matching OptCuts executable, and creates Start Menu shortcuts.

## Install

1. Obtain a setup executable from a published release or a local build.
2. Open the downloaded file. The current build is unsigned. If Windows
   SmartScreen blocks it, confirm that it came from the official release page,
   select **More info**, then select **Run anyway**.
3. Keep the setup and PowerShell windows open while the Python environment is
   prepared. A fresh installation commonly takes five to fifteen minutes.
4. Open **TopoPPI GUI** from the Start Menu.

A fresh installation needs access to GitHub, conda-forge, and
PyPI. The bundled OptCuts executable is installed locally. Once setup finishes,
routine analysis of local structures works without a network connection.

The GUI shortcut uses `pythonw.exe`, so regular desktop launches stay free of a
console window. **TopoPPI Command Prompt** opens a persistent terminal with the
dedicated environment and OptCuts ready to use. Run `topoppi --help` there to see
the commands. `TopoPPI CLI.cmd` remains in the installation folder as a direct
wrapper for scripts and automation. Both command launchers resolve paths from
their installation folder, including folders with spaces and Unicode characters.

TopoPPI and OptCuts license files are installed under
`%LOCALAPPDATA%\TopoPPI\licenses`.

## Binding analysis in 3.0

Choose **Binding analysis...** for **Binding affinity**, **Mutation effects**, or
**Rank hotspots**. Affinity supplies the first two tasks; Hotspot supplies the
third. Install the matching plugin wheels through **Plugins → Manage Plugins...
→ Install from file...**. The plugins use the
application's Python environment and include one official CPU model per task.
Their bundled parameters support offline prediction from local structures.

The [desktop guide](../../docs/desktop.md#binding-analysis) covers the three
workflows and residue-to-variant selection. The [Affinity guide](../../plugins/affinity/docs/workflow.md)
and [Hotspot guide](../../plugins/hotspot/README.md) describe outputs, model scope,
and optional custom models. Advanced neural models use their separately
installed dependencies and declared resources.

## Upgrade

Close the GUI and any TopoPPI command prompts, download the installer for the
new release, and run it. Setup updates the existing environment in
`%LOCALAPPDATA%\TopoPPI` and refreshes the launchers. Figures, manifests,
structures, and benchmark results saved elsewhere remain in place.

## Repair an installation

Run the installer for the same release again. It rechecks the Conda environment,
reinstalls the TopoPPI package, and recopies OptCuts.

If the GUI itself cannot start, it shows the location of
`%LOCALAPPDATA%\TopoPPI\gui-startup.log`. That file contains the Python startup
error.

If setup stops:

1. Read the final PowerShell message or `installation.log` in the installation
   folder for the failed download or package step.
2. Restore the network connection and run the installer again.
3. If the environment remains incomplete, uninstall TopoPPI, remove
   `%LOCALAPPDATA%\TopoPPI`, and start a fresh installation.

Setup reports the PowerShell exit code when an installation step needs attention.

## Uninstall

Close the GUI and any TopoPPI command prompts, then use **Settings > Apps >
Installed apps > TopoPPI > Uninstall** or the **Uninstall TopoPPI** Start Menu
shortcut. This removes the dedicated environment, OptCuts, launchers, and
installer files.

TopoPPI leaves analysis output in the folders you selected. The recent-file list
is stored at `%USERPROFILE%\.topoppi\gui_recent.json`; delete the `.topoppi`
folder when you also want to clear that preference data.

## Build locally

Install these tools on a 64-bit Windows development machine:

- Python 3.10, available as `python` on PATH
- Git
- CMake
- Visual Studio 2022 with the **Desktop development with C++** workload
- Inno Setup 6 or later
- network access for the pinned OptCuts source checkout

From the repository root, build OptCuts and an installer for testing this
checkout on the same computer:

```powershell
$packageSource = (Resolve-Path .).Path
$topoppiVersion = python -c "import runpy; print(runpy.run_path('src/topoppi/_version.py')['__version__'])"
.\tools\OptCuts\build_windows_optcuts.ps1 -OutputDir installer\windows
Push-Location installer\windows
& "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe" `
  "/DMyAppVersion=$topoppiVersion" `
  "/DMyPackageSpec=$packageSource" `
  TopoPPI.iss
Pop-Location
```

This local installer reads the source checkout during setup. A distributed
installer uses `MyPackageSpec` for the source archive of its matching published
commit or tag. Keep the installer version and installed package version aligned.
For the v2.1 paper release, use the build instructions stored at the v2.1 tag.

The setup executable is written to `installer/windows/Output/`.

## Packaged components

The setup executable includes `OptCuts_bin-windows-x86_64.exe`. A local build that omits it downloads the matching GitHub release asset and checks its supplied SHA-256 sidecar during setup.

The `Windows Installer` workflow builds the native executable and installer, then exercises the installed command prompt, mapping, atlas editing, export and removal. Tagged installers use that tag’s source archive. The installer and native executable are published together after the package and platform jobs succeed.
