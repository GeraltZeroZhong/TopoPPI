import hashlib
import json
import re
import unittest
from pathlib import Path

from topoppi import __version__
from topoppi.install_optcuts import LINUX_X86_64_SHA256

ROOT = Path(__file__).parents[1]


class ReleaseMetadataTests(unittest.TestCase):
    def test_public_release_metadata_uses_the_package_version(self):
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        citation_version = re.search(r'^version: "([^"]+)"$', citation, re.MULTILINE)
        citation_date = re.search(r'^date-released: "([^"]+)"$', citation, re.MULTILINE)
        self.assertIsNotNone(citation_version)
        self.assertEqual(citation_version.group(1), __version__)
        self.assertIn("cff-version: 1.2.0", citation)

        changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(f"TopoPPI {__version__}", readme)
        if f"## [{__version__}] - In preparation" in changelog:
            self.assertIsNone(citation_date, "A working-tree version has no release date yet.")
            self.assertIn("in preparation", readme)
            self.assertIn("https://github.com/GeraltZeroZhong/TopoPPI/releases", readme)
        else:
            self.assertIsNotNone(citation_date)
            self.assertIn(f"## [{__version__}] - {citation_date.group(1)}", changelog)
            self.assertIn(f"/tag/v{__version__}", readme)
        self.assertNotRegex(readme, r'(?:\]\(|src=")\./')
        self.assertIn("docs/assets/topoppi-gui-3d.png", readme)

        schema = (ROOT / "docs" / "benchmark_schema.md").read_text(encoding="utf-8")
        self.assertTrue(schema.startswith("# TopoPPI benchmark evidence schema (v2.0)\n"))
        self.assertIn("| `schema_version` | Machine-readable report layout, currently `2.0` |", schema)

        installer = (ROOT / "installer" / "windows" / "TopoPPI.iss").read_text(encoding="utf-8")
        self.assertIn(f'#define MyAppVersion "{__version__}"', installer)

        bootstrap = (ROOT / "installer" / "windows" / "install_topoppi.ps1").read_text(encoding="utf-8")
        self.assertIn(f'[string]$Version = "{__version__}"', bootstrap)
        self.assertIn("import MDAnalysis, prolif, rdkit", bootstrap)
        self.assertIn("Verifying ProLIF interaction stack", bootstrap)
        self.assertIn("micromamba-releases/releases/latest/download/micromamba-win-64", bootstrap)
        self.assertIn("Invoke-WebRequest -Uri $MicromambaUrl -OutFile $Micromamba", bootstrap)
        self.assertNotIn("tar.exe", bootstrap)
        self.assertIn('Remove-Item (Join-Path $InstallDir "TopoPPI GUI.cmd")', bootstrap)

        for platform_name in ("windows", "macos"):
            platform_readme = (ROOT / "installer" / platform_name / "README.md").read_text(encoding="utf-8")
            self.assertIn(f"TopoPPI {__version__}", platform_readme)

    def test_windows_installer_propagates_bootstrap_failure_and_hides_the_gui_console(self):
        installer = (ROOT / "installer" / "windows" / "TopoPPI.iss").read_text(encoding="utf-8")
        launcher = (ROOT / "installer" / "windows" / "launch_gui.pyw").read_text(encoding="utf-8")
        uninstaller = (ROOT / "installer" / "windows" / "uninstall_topoppi.ps1").read_text(encoding="utf-8")
        manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")

        self.assertIn("CurStepChanged", installer)
        self.assertIn("ResultCode <> 0", installer)
        self.assertIn("RaiseException", installer)
        self.assertIn("GetCustomSetupExitCode", installer)
        self.assertIn("BootstrapExitCode := ResultCode", installer)
        self.assertIn(r'Filename: "{app}\env\pythonw.exe"', installer)
        self.assertIn("launch_gui.pyw", installer)
        self.assertNotIn("TopoPPI GUI.cmd", installer)
        self.assertIn(r'Name: "{group}\TopoPPI Command Prompt"', installer)
        self.assertIn(r'Filename: "{app}\TopoPPI Command Prompt.cmd"', installer)
        self.assertIn(r'Name: "{group}\TopoPPI CLI.lnk"', installer)
        self.assertIn('os.environ["TOPOPPI_OPTCUTS_BIN"]', launcher)
        self.assertIn("TopoPPI GUI.cmd", uninstaller)
        self.assertIn("TopoPPI CLI.cmd", uninstaller)
        self.assertIn("TopoPPI Command Prompt.cmd", uninstaller)
        self.assertIn("gui-startup.log", uninstaller)
        self.assertIn("installation.log", uninstaller)
        self.assertIn('Type: dirifempty; Name: "{app}"', installer)
        self.assertIn("TopoPPI-LICENSE.txt", installer)
        self.assertIn("OptCuts-LICENSE.txt", installer)
        self.assertIn("OptCuts-NOTICE.md", installer)
        self.assertIn("OptCuts-THIRD-PARTY-LICENSES.txt", installer)
        self.assertIn("*.pyw", manifest)
        self.assertIn("docs/assets/topoppi-gui-sanitized.png", manifest)
        self.assertNotIn("3ff22687bd403a67cd66caeacc95baee.png", manifest)

    @unittest.skipUnless((ROOT / ".gitignore").is_file(), ".gitignore is not included in source distributions")
    def test_generated_prolif_manifests_are_ignored(self):
        gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")

        self.assertIn("*.prolif.csv", gitignore.splitlines())

        bootstrap = (ROOT / "installer" / "windows" / "install_topoppi.ps1").read_text(encoding="utf-8")
        self.assertIn('Join-Path $InstallDir "TopoPPI CLI.cmd"', bootstrap)
        self.assertIn('Join-Path $InstallDir "TopoPPI Command Prompt.cmd"', bootstrap)
        self.assertIn('set "PATH=%~dp0env\\Scripts;%~dp0env;%PATH%"', bootstrap)
        self.assertIn('"%~dp0env\\Scripts\\topoppi.exe" %*', bootstrap)
        self.assertIn("cmd.exe /K", bootstrap)

    def test_macos_launcher_explains_first_start_and_recovery(self):
        launcher = (ROOT / "installer" / "macos" / "TopoPPI").read_text(encoding="utf-8")

        self.assertIn("show_startup_error", launcher)
        self.assertIn("display dialog", launcher)
        self.assertIn("launcher.log", launcher)
        self.assertIn('rm -rf "$environment_dir"', launcher)

        build_script = (ROOT / "installer" / "macos" / "build_app.sh").read_text(encoding="utf-8")
        self.assertIn("TopoPPI-LICENSE.txt", build_script)
        self.assertIn("OptCuts-LICENSE.txt", build_script)
        self.assertIn("OptCuts-NOTICE.md", build_script)
        self.assertIn("OptCuts-THIRD-PARTY-LICENSES.txt", build_script)
        self.assertIn(
            "--exclude 'lib/python*/site-packages/topoppi-*.dist-info/direct_url.json'",
            build_script,
        )
        self.assertNotIn("direct_url.json' -delete", build_script)
        self.assertNotIn('find "$environment_prefix" -type f -path', build_script)

    def test_formal_example_uses_the_release_optcuts_artifact(self):
        payload = json.loads((ROOT / "docs" / "benchmark_config.example.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["optcuts"]["expected_binary_sha256"], LINUX_X86_64_SHA256)
        self.assertGreater(payload["optcuts"]["residue_fragmentation_weight"], 0.0)
        self.assertIn("residue_aware_optcuts", payload["optcuts_variants"])

    @unittest.skipUnless(
        (ROOT / ".github" / "workflows" / "publish.yml").is_file(),
        "GitHub workflows are not included in source distributions",
    )
    def test_release_jobs_bind_gh_to_the_target_repository(self):
        workflow = (ROOT / ".github" / "workflows" / "publish.yml").read_text(encoding="utf-8")

        self.assertEqual(workflow.count("GH_REPO: ${{ github.repository }}"), 2)

    @unittest.skipUnless(
        (ROOT / ".github" / "workflows" / "windows-installer.yml").is_file(),
        "GitHub workflows are not included in source distributions",
    )
    def test_native_installer_workflows_construct_the_desktop_app(self):
        windows = (ROOT / ".github" / "workflows" / "windows-installer.yml").read_text(encoding="utf-8")
        macos = (ROOT / ".github" / "workflows" / "macos-app.yml").read_text(encoding="utf-8")

        for workflow in (windows, macos):
            self.assertIn("from topoppi.gui_app import ProtSurfApp", workflow)
            self.assertIn("ProtSurfApp(root)", workflow)
            self.assertIn("root.update_idletasks()", workflow)
        self.assertIn("TopoPPI Command Prompt.lnk", windows)
        self.assertIn("& $cli --version", windows)
        self.assertIn("Start-Process -FilePath $installer", windows)
        self.assertIn("Start-Process -FilePath $uninstaller", windows)
        self.assertEqual(windows.count("| Wait-Process -Timeout"), 1)
        self.assertEqual(windows.count(") -Wait -PassThru"), 1)
        self.assertIn("timeout-minutes: 40", windows)
        self.assertIn('"/DIR=`"$installDir`""', windows)
        self.assertIn('"/LOG=`"$installLog`""', windows)
        self.assertIn('"TopoPPI 用户"', windows)
        self.assertIn('@("topoppi --version", "exit") | & $commandPrompt', windows)
        self.assertIn(r"tests\fixtures\tiny_complex.pdb", windows)
        self.assertNotIn(r"tests\fixtures\1bvk.pdb", windows)
        self.assertIn("tests/fixtures/1bvk.pdb", macos)
        self.assertIn('"cmake<4"', macos)
        self.assertEqual(macos.count('export PATH="$RUNNER_TEMP/topoppi-macos-build/bin:$PATH"'), 2)
        for workflow in (windows, macos):
            self.assertIn("--map-style footprints --export-atlas", workflow)
            self.assertIn("smoke_footprints.py", workflow)

    @unittest.skipUnless(
        (ROOT / "tools" / "OptCuts" / "OptCuts_bin").is_file(),
        "OptCuts binary is not included in source distributions",
    )
    def test_tracked_linux_binary_matches_the_release_digest(self):
        binary = ROOT / "tools" / "OptCuts" / "OptCuts_bin"

        self.assertEqual(hashlib.sha256(binary.read_bytes()).hexdigest(), LINUX_X86_64_SHA256)


if __name__ == "__main__":
    unittest.main()
