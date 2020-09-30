
PUSHD "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build"
::PUSHD F:\VS19\VS_IDE\VC\Auxiliary\Build
call vcvars64.bat
POPD

::curl -O https://dl.bintray.com/boostorg/release/1.74.0/source/boost_1_74_0.zip

PUSHD ".\boost_1_74_0\boost_1_74_0"
call bootstrap.bat
::copy /Y ..\..\files\project-config.jam .
.\b2 headers
.\b2 -q -a -sZLIB_SOURCE="..\..\zlib-1.2.11" -sZLIB_INCLUDE="..\..\zlib-1.2.11" -sZLIB_LIBPATH="..\..\zlib-1.2.11" --toolset=msvc-14.2 address-model=64 runtime-link=static threading=multi stage
::.\b2 -q -a --with-iostreams --toolset=msvc-14.2 address-model=64 runtime-link=shared threading=multi stage
::@REM bjam --toolset=msvc-12.0 address-model=64 runtime-link=static runtime-debugging=on variant=debug threading=multi stage
::.\b2 -q -a --toolset=msvc-14.2 address-model=64 runtime-link=static threading=single stage
::@REM bjam --toolset=msvc-12.0 address-model=64 runtime-link=static runtime-debugging=on variant=debug threading=single stage
POPD
