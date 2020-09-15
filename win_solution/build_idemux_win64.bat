PUSHD "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build"
::PUSHD F:\VS19\VS_IDE\VC\Auxiliary\Build
call vcvars64.bat
POPD

::.\build_zlib.bat

::.\build_boost.bat

::"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\MSBuild\Current\Bin\amd64\msbuild.exe" "idemuxcpp.sln" /t:Rebuild /p:configuration=Release-cluster

if not exist "iDemux_win10_64bit" mkdir iDemux_win10_64bit
if not exist "iDemux_win10_64bit\tests" mkdir iDemux_win10_64bit\tests
if not exist "iDemux_win10_64bit\bin" mkdir iDemux_win10_64bit\bin
if not exist "iDemux_win10_64bit\share" mkdir iDemux_win10_64bit\share
if not exist "iDemux_win10_64bit\share\idemuxCPP" mkdir iDemux_win10_64bit\share\idemuxcpp
Xcopy .\x64\Release-cluster\idemuxcpp.exe iDemux_win10_64bit\bin\
Xcopy .\x64\Release-cluster\test_parser.exe iDemux_win10_64bit\tests\
Xcopy .\x64\Release-cluster\test_demuxer.exe iDemux_win10_64bit\tests\
Xcopy .\x64\Release-cluster\test_barcode.exe iDemux_win10_64bit\tests\
Xcopy ..\LICENCE iDemux_win10_64bit\LICENCE*.txt
Xcopy ..\README.md iDemux_win10_64bit\

:: copy boost support dll or let user install the package Visual C++ Redistributable for Visual Studio 2015 (https://www.microsoft.com/de-at/download/details.aspx?id=48145).
Xcopy .\vcomp140.dll iDemux_win10_64bit\tests\
Xcopy .\vcomp140.dll iDemux_win10_64bit\bin\
if not exist "iDemux_win10_64bit\tests\resources" Xcopy ..\tests\resources .\iDemux_win10_64bit\tests\resources /E/H/C/I
if not exist "iDemux_win10_64bit\share\idemuxCPP\barcodes" Xcopy ..\misc\barcodes iDemux_win10_64bit\share\idemuxcpp\barcodes /E/H/C/I
