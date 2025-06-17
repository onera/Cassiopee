; -- CassiopeeWin32.iss --
; Inno Setup file

[Setup]
AppName=Cassiopee
AppVersion=4.1
DefaultDirName={code:DefDirRoot}\Cassiopee
DefaultGroupName=Cassiopee
Compression=lzma2
SolidCompression=yes
OutputBaseFilename=Cassiopee-4.1-win32
PrivilegesRequired=lowest
AppPublisher=ONERA

[Files]
Source: "D:\Documents and Settings\gleize\Mes documents\Cassiopee\Dist\*"; DestDir: "{app}\Dist" ; Flags: recursesubdirs createallsubdirs
Source: "D:\Documents and Settings\gleize\Mes documents\Cassiopee\env_Cassiopee_win32.bat"; DestDir: "{app}\Dist"; AfterInstall: CurStepChanged()

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
  fileName := ExpandConstant('{app}\Dist\env_Cassiopee_win32.bat');
  CassiopeeVar := ExpandConstant('{app}');
  SetArrayLength(lines, 6);
  lines[0] := 'set CASSIOPEE='+CassiopeeVar;
  lines[1] := 'path = %CASSIOPEE%\Dist\bin\x86\Python27;%CASSIOPEE%\Dist\bin\x86\Python27\Scripts;%CASSIOPEE%\Dist\bin\x86\Python27\DLLs;%PATH%';
  lines[2] := 'path = %CASSIOPEE%\Dist\bin\x86;%CASSIOPEE%\Dist\bin\x86\Lib;%CASSIOPEE%\Dist\bin\x86\Lib\site-packages;%PATH%';
  lines[3] := 'set PYTHONHOME=%CASSIOPEE%\Dist\bin\x86\Python27';
  lines[4] := 'set ELSAPROD=x86';
  lines[5] := 'set PYTHONPATH=%CASSIOPEE%\Dist\bin\x86;%CASSIOPEE%\Dist\bin\x86\Python27\Scripts;%CASSIOPEE%\Dist\bin\x86\Lib\site-packages';
  SaveStringsToFile(filename,lines,false);

  fileName := ExpandConstant('{app}\Dist\bin\x86\cassiopeeRunWin32.bat');
  CassiopeeVar := ExpandConstant('{app}');
  SetArrayLength(lines, 7);
  lines[0] := 'set CASSIOPEE='+CassiopeeVar;
  lines[1] := 'path = %CASSIOPEE%\Dist\bin\x86\Python27;%CASSIOPEE%\Dist\bin\x86\Python27\Scripts;%CASSIOPEE%\Dist\bin\x86\Python27\DLLs;%PATH%';
  lines[2] := 'path = %CASSIOPEE%\Dist\bin\x86;%CASSIOPEE%\Dist\bin\x86\Lib;%CASSIOPEE%\Dist\bin\x86\Lib\site-packages;%PATH%';
  lines[3] := 'set PYTHONHOME=%CASSIOPEE%\Dist\bin\x86\Python27';
  lines[4] := 'set ELSAPROD=x86';
  lines[5] := 'set PYTHONPATH=%CASSIOPEE%\Dist\bin\x86;%CASSIOPEE%\Dist\bin\x86\Python27\Scripts;%CASSIOPEE%\Dist\bin\x86\Lib\site-packages'
  lines[6] := 'python "%CASSIOPEE%\Dist\bin\x86\tkCassiopee.pyc" %1'
  Result := SaveStringsToFile(filename,lines,false);
  exit;
end;
 
procedure CurStepChanged();
begin
    CreateBatch();
end;

[Icons]
Name: "{group}\Cassiopee"; Filename: "{app}\Dist\bin\x86\cassiopeeRunWin32.bat" ; Flags: runminimized ; WorkingDir: "%USERPROFILE%" ; IconFilename: "{app}\Dist\bin\x86\Lib\site-packages\CPlot\logoCassiopee32.ico"
Name: "{group}\Command shell"; Filename: "cmd.exe" ; Parameters: "/k ""{app}\Dist\env_Cassiopee_win32.bat""" ; WorkingDir: "%USERPROFILE%" ; Flags: runmaximized
Name: "{group}\Uninstall"; Filename: "{uninstallexe}"
