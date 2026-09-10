#define MyAppName "TopoPPI"
#ifndef MyAppVersion
#define MyAppVersion "2.1"
#endif
#ifndef MyPackageSpec
#define MyPackageSpec ""
#endif

[Setup]
AppId={{4F4E6672-FF5E-43AF-9DC5-5E3E64CE4FE3}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher=TopoPPI contributors
AppPublisherURL=https://github.com/GeraltZeroZhong/TopoPPI
AppSupportURL=https://github.com/GeraltZeroZhong/TopoPPI/issues
AppUpdatesURL=https://github.com/GeraltZeroZhong/TopoPPI/releases
DefaultDirName={localappdata}\TopoPPI
DefaultGroupName=TopoPPI
DisableProgramGroupPage=yes
OutputDir=Output
OutputBaseFilename=TopoPPI-{#MyAppVersion}-windows-x86_64-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
SetupIconFile=..\..\src\topoppi\assets\topoppi.ico
UninstallDisplayIcon={app}\installer\assets\topoppi.ico

[Files]
Source: "install_topoppi.ps1"; DestDir: "{app}\installer"; Flags: ignoreversion
Source: "uninstall_topoppi.ps1"; DestDir: "{app}\installer"; Flags: ignoreversion
Source: "launch_gui.pyw"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\src\topoppi\assets\topoppi.ico"; DestDir: "{app}\installer\assets"; Flags: ignoreversion
Source: "..\..\LICENSE"; DestDir: "{app}\licenses"; DestName: "TopoPPI-LICENSE.txt"; Flags: ignoreversion
Source: "..\..\tools\OptCuts\LICENSE.txt"; DestDir: "{app}\licenses"; DestName: "OptCuts-LICENSE.txt"; Flags: ignoreversion
Source: "..\..\tools\OptCuts\NOTICE.md"; DestDir: "{app}\licenses"; DestName: "OptCuts-NOTICE.md"; Flags: ignoreversion
Source: "..\..\tools\OptCuts\THIRD_PARTY_LICENSES.txt"; DestDir: "{app}\licenses"; DestName: "OptCuts-THIRD-PARTY-LICENSES.txt"; Flags: ignoreversion
#ifexist "OptCuts_bin-windows-x86_64.exe"
Source: "OptCuts_bin-windows-x86_64.exe"; DestDir: "{app}\installer\assets"; Flags: ignoreversion
#endif

[Icons]
Name: "{group}\TopoPPI GUI"; Filename: "{app}\env\pythonw.exe"; Parameters: """{app}\launch_gui.pyw"""; WorkingDir: "{app}"; IconFilename: "{app}\installer\assets\topoppi.ico"
Name: "{group}\TopoPPI Command Prompt"; Filename: "{app}\TopoPPI Command Prompt.cmd"; WorkingDir: "{app}"; IconFilename: "{app}\installer\assets\topoppi.ico"
Name: "{group}\Uninstall TopoPPI"; Filename: "{uninstallexe}"; IconFilename: "{app}\installer\assets\topoppi.ico"

[InstallDelete]
Type: files; Name: "{group}\TopoPPI CLI.lnk"

[UninstallRun]
Filename: "powershell.exe"; Parameters: "-NoProfile -ExecutionPolicy Bypass -File ""{app}\installer\uninstall_topoppi.ps1"" -InstallDir ""{app}"""; Flags: waituntilterminated runhidden

[UninstallDelete]
Type: dirifempty; Name: "{app}"

[Code]
var
  BootstrapExitCode: Integer;

function GetCustomSetupExitCode: Integer;
begin
  Result := BootstrapExitCode;
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  Parameters: String;
  ResultCode: Integer;
begin
  if CurStep <> ssPostInstall then
    Exit;

  Parameters :=
    '-NoProfile -ExecutionPolicy Bypass -File "' +
    ExpandConstant('{app}\installer\install_topoppi.ps1') +
    '" -InstallDir "' + ExpandConstant('{app}') +
    '" -Version "{#MyAppVersion}" -PackageSpec "{#MyPackageSpec}"';

  if not Exec(
    ExpandConstant('{sys}\WindowsPowerShell\v1.0\powershell.exe'),
    Parameters,
    ExpandConstant('{app}'),
    SW_SHOWNORMAL,
    ewWaitUntilTerminated,
    ResultCode
  ) then
  begin
    BootstrapExitCode := 1;
    RaiseException('TopoPPI could not start its environment setup.');
  end;

  if ResultCode <> 0 then
  begin
    BootstrapExitCode := ResultCode;
    RaiseException(
      'TopoPPI setup failed with exit code ' + IntToStr(ResultCode) +
      '. See ' + ExpandConstant('{app}\installation.log') + ' for details.'
    );
  end;
end;
