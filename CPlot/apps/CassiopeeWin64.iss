; -- CassiopeeWin64.iss --
; Inno Setup file

[Setup]
AppName=Cassiopee
AppVersion=3.3
DefaultDirName={code:DefDirRoot}\Cassiopee
DefaultGroupName=Cassiopee
Compression=lzma2
SolidCompression=yes
OutputBaseFilename=Cassiopee-3.3-win64
PrivilegesRequired=lowest
AppPublisher=ONERA

[Files]
Source: "D:\benoit\Cassiopee\Dist\*"; DestDir: "{app}\Dist" ; Flags: recursesubdirs createallsubdirs
Source: "D:\benoit\Cassiopee\env_Cassiopee_win64.bat"; DestDir: "{app}\Dist"; AfterInstall: CurStepChanged()

[Code]
function IsRegularUser(): boolean;
begin
Result := not (isAdminLoggedOn or IsPowerUserLoggedOn);
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
  lines[1] := 'path=%CASSIOPEE%\Dist\bin\win64\lib;%CASSIOPEE%\Dist\bin\win64\bin;%CASSIOPEE%\Dist\bin\win64;%PATH%';
  lines[2] := 'set PYTHONEXE=python3.8';
  lines[3] := 'set PYTHONHOME=%CASSIOPEE%\Dist\bin\win64';
  lines[4] := 'set ELSAPROD=win64';
  lines[5] := 'set PYTHONPATH=%CASSIOPEE%\Dist\bin\win64\bin;%CASSIOPEE%\Dist\bin\win64\lib\python3.8;%CASSIOPEE%\Dist\bin\win64\lib\python3.8\site-packages';
  lines[6] := 'set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%'
  SaveStringsToFile(filename,lines,false);

  fileName := ExpandConstant('{app}\Dist\bin\win64\cassiopeeRunWin64.bat');
  CassiopeeVar := ExpandConstant('{app}');
  SetArrayLength(lines, 8);
  lines[0] := 'set CASSIOPEE='+CassiopeeVar;
  lines[1] := 'path=%CASSIOPEE%\Dist\bin\win64\bin;%CASSIOPEE%\Dist\bin\win64;%PATH%';
  lines[2] := 'set PYTHONEXE=python3.8';
  lines[3] := 'set PYTHONHOME=%CASSIOPEE%\Dist\bin\win64';
  lines[4] := 'set ELSAPROD=win64';
  lines[5] := 'set PYTHONPATH=%CASSIOPEE%\Dist\bin\win64;%CASSIOPEE%\Dist\bin\win64\Lib\python3.8;%CASSIOPEE%\Dist\bin\win64\lib\python3.8\site-packages'
  lines[6] := 'python "%CASSIOPEE%\Dist\bin\win64\tkCassiopee.py" %1'
  lines[7] := 'set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%'
  Result := SaveStringsToFile(filename,lines,false);
  exit;
end;
 
procedure CurStepChanged();
begin
    CreateBatch();
end;

[Icons]
Name: "{group}\Cassiopee"; Filename: "{app}\Dist\bin\win64\cassiopeeRunWin64.bat" ; Flags: runminimized ; WorkingDir: "%USERPROFILE%" ; IconFilename: "{app}\Dist\bin\win64\lib\python3.8\site-packages\CPlot\logoCassiopee32.ico"
Name: "{group}\Command shell"; Filename: "cmd.exe" ; Parameters: "/k ""{app}\Dist\env_Cassiopee_win64.bat""" ; WorkingDir: "%USERPROFILE%" ; Flags: runmaximized
Name: "{group}\Uninstall"; Filename: "{uninstallexe}"
