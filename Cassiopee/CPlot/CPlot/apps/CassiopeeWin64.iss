; -- CassiopeeWin64.iss --
; Inno Setup file

[Setup]
AppName=Cassiopee
AppVersion=2026.11
DefaultDirName={code:DefDirRoot}\Cassiopee
DefaultGroupName=Cassiopee
Compression=lzma2
SolidCompression=yes
OutputBaseFilename=Cassiopee-2026.11-win64
PrivilegesRequired=lowest
AppPublisher=ONERA

[Files]
Source: "C:\Users\benoit\Code\Cassiopee\Dist\*"; DestDir: "{app}\Dist" ; Flags: recursesubdirs createallsubdirs; Excludes: "ubuntu"
Source: "C:\Users\benoit\Code\Cassiopee\Cassiopee\Envs\env_Cassiopee_win64.bat"; DestDir: "{app}\Dist"; AfterInstall: CurStepChanged()
Source: "C:\Users\benoit\Code\Cassiopee\Cassiopee\CPlot\CPlot\apps\cassiopeeRunWin64.bat"; DestDir: "{app}\Dist"; AfterInstall: CurStepChanged()

[Code]
function IsRegularUser(): boolean;
begin
Result := not (isAdmin or IsPowerUserLoggedOn);
end;
function DefDirRoot(Param: String): String;
begin
if IsRegularUser then
Result := ExpandConstant('{localappdata}')
else
Result := ExpandConstant('{pf}')
end;
function CreateBatch(): boolean;
var
  fileName : string;
  CassiopeeVar : string;
  lines : TArrayOfString;
begin
  Result := true;
  fileName := ExpandConstant('{app}\Dist\env_Cassiopee_win64.bat');
  CassiopeeVar := ExpandConstant('{app}');
  SetArrayLength(lines, 7);
  lines[0] := 'set CASSIOPEE='+CassiopeeVar;
  lines[1] := 'path=%CASSIOPEE%\Dist\bin\win64;%CASSIOPEE%\Dist\bin\win64\lib;%CASSIOPEE%\Dist\bin\win64\lib2;%CASSIOPEE%\Dist\bin\win64\bin2;%CASSIOPEE%\Dist\bin\win64\bin;%PATH%';
  lines[2] := 'set PYTHONEXE=python3.12';
  lines[3] := 'set PYTHONHOME=%CASSIOPEE%\Dist\bin\win64';
  lines[4] := 'set ELSAPROD=win64';
  lines[5] := 'set PYTHONPATH=%CASSIOPEE%\Dist\bin\win64\lib\python3.12;%CASSIOPEE%\Dist\bin\win64\lib\python3.12\site-packages;%CASSIOPEE%\Dist\bin\win64\lib2\python3.12;%CASSIOPEE%\Dist\bin\win64\lib2\python3.12\site-packages;%CASSIOPEE%\Dist\bin\win64\lib2\python3.12\lib-dynload';
  lines[6] := 'set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%';
  SaveStringsToFile(filename,lines,false);

  fileName := ExpandConstant('{app}\Dist\env_Cassiopee_win64.ps1');
  CassiopeeVar := ExpandConstant('{app}');
  SetArrayLength(lines, 7);
  lines[0] := '$env:CASSIOPEE = "'+CassiopeeVar+'"';
  lines[1] := '$env:PATH = "$env:CASSIOPEE\Dist\bin\win64;$env:CASSIOPEE\Dist\bin\win64\lib;$env:CASSIOPEE\Dist\bin\win64\lib2;$env:CASSIOPEE\Dist\bin\win64\bin2;$env:CASSIOPEE\Dist\bin\win64\bin;$env:PATH"';
  lines[2] := '$env:PYTHONEXE = "python3.12"';
  lines[3] := '$env:PYTHONHOME = "$env:CASSIOPEE\Dist\bin\win64"';
  lines[4] := '$env:ELSAPROD = "win64"';
  lines[5] := '$env:PYTHONPATH = "$env:CASSIOPEE\Dist\bin\win64\lib\python3.12;$env:CASSIOPEE\Dist\bin\win64\lib\python3.12\site-packages;$env:CASSIOPEE\Dist\bin\win64\lib2\python3.12;$env:CASSIOPEE\Dist\bin\win64\lib2\python3.12\site-packages;$env:CASSIOPEE\Dist\bin\win64\lib2\python3.12\lib-dynload"'
  lines[6] := '$env:OMP_NUM_THREADS = $env:NUMBER_OF_PROCESSORS';
  SaveStringsToFile(filename,lines,false);

  fileName := ExpandConstant('{app}\Dist\cassiopeeRunWin64.bat');
  CassiopeeVar := ExpandConstant('{app}');
  SetArrayLength(lines, 8);
  lines[0] := 'set CASSIOPEE='+CassiopeeVar;
  lines[1] := 'path=%CASSIOPEE%\Dist\bin\win64;%CASSIOPEE%\Dist\bin\win64\lib;%CASSIOPEE%\Dist\bin\win64\lib2;%CASSIOPEE%\Dist\bin\win64\bin2;%CASSIOPEE%\Dist\bin\win64\bin;%PATH%';
  lines[2] := 'set PYTHONEXE=python3.12';
  lines[3] := 'set PYTHONHOME=%CASSIOPEE%\Dist\bin\win64';
  lines[4] := 'set ELSAPROD=win64';
  lines[5] := 'set PYTHONPATH=%CASSIOPEE%\Dist\bin\win64\lib\python3.12;%CASSIOPEE%\Dist\bin\win64\lib\python3.12\site-packages;%CASSIOPEE%\Dist\bin\win64\lib2\python3.12;%CASSIOPEE%\Dist\bin\win64\lib2\python3.12\site-packages;%CASSIOPEE%\Dist\bin\win64\lib2\python3.12\lib-dynload';
  lines[6] := 'python "%CASSIOPEE%\Dist\bin\win64\tkCassiopee.py" %1';
  lines[7] := 'set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%';
  Result := SaveStringsToFile(filename,lines,false);
  
  exit;
end;
 
procedure CurStepChanged();
begin
    CreateBatch();
end;

[Icons]
Name: "{group}\Cassiopee"; Filename: "cmd.exe" ; Parameters: "/k ""{app}\Dist\cassiopeeRunWin64.bat""" ; Flags: runminimized ; WorkingDir: "%USERPROFILE%" ; IconFilename: "{app}\Dist\bin\win64\lib\python3.12\site-packages\CPlot\logoCassiopee32.ico"
Name: "{group}\Dos shell"; Filename: "cmd.exe" ; Parameters: "/k ""{app}\Dist\env_Cassiopee_win64.bat""" ; WorkingDir: "%USERPROFILE%" ; Flags: runmaximized
Name: "{group}\Power shell"; Filename: "powershell.exe" ; Parameters: "-ExecutionPolicy Bypass -noexit ""{app}\Dist\env_Cassiopee_win64.ps1""" ; WorkingDir: "%USERPROFILE%" ; Flags: runmaximized
Name: "{group}\Uninstall"; Filename: "{uninstallexe}"
