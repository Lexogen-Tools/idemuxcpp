
::PUSHD "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build"
PUSHD F:\VS19\VS_IDE\VC\Auxiliary\Build
call vcvars64.bat
POPD

PUSHD ".\boost"
call bootstrap.bat
::copy /Y ..\..\files\project-config.jam .
.\b2 headers
::bjam --toolset=msvc-14.2 address-model=64 runtime-link=static threading=multi stage
::@REM bjam --toolset=msvc-12.0 address-model=64 runtime-link=static runtime-debugging=on variant=debug threading=multi stage
.\b2 --toolset=msvc-14.2 address-model=64 runtime-link=static threading=single stage
::@REM bjam --toolset=msvc-12.0 address-model=64 runtime-link=static runtime-debugging=on variant=debug threading=single stage
POPD
