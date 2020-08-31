:: build tools for visual c++: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019

::"C:\Program Files (x86)\Microsoft Visual Studio 14.0\"

::PUSHD "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build"
PUSHD F:\VS19\VS_IDE\VC\Auxiliary\Build
call vcvars64.bat
POPD

curl -O https://zlib.net/zlib-1.2.11.tar.gz
tar -xzf zlib-1.2.11.tar.gz
::TODO: change cl compiler flag from /MD to /MT (multithreads and static lib)
pushd zlib-1.2.11
::"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\"
nmake -f win32/Makefile.msc
popd
