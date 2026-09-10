[CmdletBinding()]
param(
    [string]$InstallDir = "$env:LOCALAPPDATA\TopoPPI",
    [string]$Version = "2.1",
    [string]$PackageSpec = "",
    [string]$MicromambaUrl = "https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-win-64"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest
$TranscriptStarted = $false
$CacheDrive = $null
$Subst = Join-Path $env:SystemRoot "System32\subst.exe"

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host "==> $Message"
}

function Invoke-External {
    param(
        [Parameter(Mandatory = $true)][string]$FilePath,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )
    & $FilePath @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$FilePath failed with exit code $LASTEXITCODE"
    }
}

try {
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

    $InstallDir = [System.IO.Path]::GetFullPath($InstallDir)
    $BinDir = Join-Path $InstallDir "bin"
    $EnvDir = Join-Path $InstallDir "env"
    $RootPrefix = Join-Path $InstallDir "mamba-root"
    $OptCutsDir = Join-Path $InstallDir "bin"
    $Micromamba = Join-Path $BinDir "micromamba.exe"
    $BundledOptCuts = Join-Path $InstallDir "installer\assets\OptCuts_bin-windows-x86_64.exe"

    New-Item -ItemType Directory -Force -Path $InstallDir, $BinDir | Out-Null
    Start-Transcript -Path (Join-Path $InstallDir "installation.log") -Force | Out-Null
    $TranscriptStarted = $true

    if (!(Test-Path $Micromamba)) {
        Write-Step "Downloading micromamba"
        Invoke-WebRequest -Uri $MicromambaUrl -OutFile $Micromamba
    }

    # Micromamba's index reader and archive extractor need an ASCII cache path.
    # Keep the environment at its selected path and map only cache access during setup.
    if ($RootPrefix -match '[^\x00-\x7F]') {
        $UsedDrives = [System.IO.Directory]::GetLogicalDrives()
        foreach ($Code in 90..68) {
            $CandidateDrive = "{0}:" -f [char]$Code
            if ($UsedDrives -notcontains "$CandidateDrive\") {
                Invoke-External $Subst @($CandidateDrive, $InstallDir)
                $CacheDrive = $CandidateDrive
                $RootPrefix = "$CacheDrive\$([System.IO.Path]::GetRandomFileName())"
                break
            }
        }
        if ($null -eq $CacheDrive) {
            throw "Windows setup needs an unused drive letter for its package cache."
        }
    }
    $env:MAMBA_ROOT_PREFIX = $RootPrefix
    $env:PYTHONUTF8 = "1"
    Invoke-External $Micromamba @("--version")
    $Packages = @(
        "python=3.10",
        "tk",
        "igl=2.6.*",
        "numpy>=1.21",
        "scipy",
        "biopython",
        "scikit-image",
        "matplotlib",
        "trimesh",
        "networkx",
        "pillow",
        "rtree",
        "shapely",
        "mdanalysis",
        "rdkit",
        "psutil",
        "tqdm",
        "meshio",
        "pip"
    )

    if (Test-Path $EnvDir) {
        Write-Step "Updating TopoPPI Conda environment"
        Invoke-External $Micromamba (@("install", "-y", "-p", $EnvDir, "-c", "conda-forge") + $Packages)
    }
    else {
        Write-Step "Creating TopoPPI Conda environment"
        Invoke-External $Micromamba (@("create", "-y", "-p", $EnvDir, "-c", "conda-forge") + $Packages)
    }

    $Python = Join-Path $EnvDir "python.exe"
    $OptCutsInstaller = Join-Path $EnvDir "Scripts\topoppi-install-optcuts.exe"
    $OptCutsExe = Join-Path $OptCutsDir "OptCuts_bin.exe"
    if (!(Test-Path $Python)) {
        throw "Python was not installed at $Python."
    }

    $ResolvedPackageSpec = $PackageSpec
    if ([string]::IsNullOrWhiteSpace($ResolvedPackageSpec)) {
        $ResolvedPackageSpec = "topoppi==$Version"
    }

    Write-Step "Installing TopoPPI from $ResolvedPackageSpec"
    Invoke-External $Python @("-m", "pip", "install", "--upgrade", "--force-reinstall", "pip")
    Invoke-External $Python @("-m", "pip", "install", "prolif>=2.0")
    Invoke-External $Python @("-m", "pip", "install", "--no-deps", "--force-reinstall", $ResolvedPackageSpec)

    Write-Step "Verifying ProLIF interaction stack"
    Invoke-External $Python @(
        "-c",
        "import MDAnalysis, prolif, rdkit; print('ProLIF interaction stack ready')"
    )

    Write-Step "Installing Windows OptCuts artifact"
    if (Test-Path $BundledOptCuts) {
        Copy-Item -Force $BundledOptCuts $OptCutsExe
    }
    else {
        if (!(Test-Path $OptCutsInstaller)) {
            throw "topoppi-install-optcuts was not installed at $OptCutsInstaller."
        }
        Invoke-External $OptCutsInstaller @("--platform", "windows-x86_64", "--install-dir", $OptCutsDir, "--force")
    }

    if (!(Test-Path $OptCutsExe)) {
        throw "OptCuts was not installed at $OptCutsExe."
    }

    Write-Step "Writing launchers"
    Remove-Item (Join-Path $InstallDir "TopoPPI GUI.cmd") -Force -ErrorAction SilentlyContinue
    $CliLauncher = @"
@echo off
set "PYTHONUTF8=1"
set "TOPOPPI_HOME=%~dp0"
set "TOPOPPI_OPTCUTS_BIN=%~dp0bin\OptCuts_bin.exe"
"%~dp0env\Scripts\topoppi.exe" %*
"@
    Set-Content -Path (Join-Path $InstallDir "TopoPPI CLI.cmd") -Value $CliLauncher -Encoding ASCII

    $CommandPromptLauncher = @"
@echo off
set "PYTHONUTF8=1"
set "TOPOPPI_HOME=%~dp0"
set "TOPOPPI_OPTCUTS_BIN=%~dp0bin\OptCuts_bin.exe"
set "PATH=%~dp0env\Scripts;%~dp0env;%PATH%"
cd /d "%~dp0"
echo TopoPPI $Version command prompt
echo Run topoppi --help to see the available commands.
echo Run exit or close this window when you are finished.
cmd.exe /K
"@
    Set-Content -Path (Join-Path $InstallDir "TopoPPI Command Prompt.cmd") -Value $CommandPromptLauncher -Encoding ASCII

    Write-Step "TopoPPI installation finished"
}
catch {
    Write-Error "TopoPPI installation failed: $($_.Exception.Message)" -ErrorAction Continue
    Write-Host "Installation details: $(Join-Path $InstallDir 'installation.log')"
    exit 1
}
finally {
    if ($null -ne $CacheDrive) {
        if (Test-Path $RootPrefix) {
            Remove-Item -Recurse -Force $RootPrefix
        }
        Invoke-External $Subst @($CacheDrive, "/D")
    }
    if ($TranscriptStarted) {
        Stop-Transcript | Out-Null
    }
}
