@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d %~dp0
set INC=/I ..\include /I D:\GoLD_rayt\GoLD_rayt\include /I D:\GoLD_rayt\GoLD_rayt\external
set GGX=D:\GoLD_rayt\GoLD_rayt\src\Microfacet\GGX.cpp
echo === RGB build
cl /nologo /O2 /std:c++20 /EHsc /W4 /utf-8 /openmp /DNDEBUG %INC% compile_check.cpp %GGX% /Fe:check_rgb.exe
echo === spectral build
cl /nologo /O2 /std:c++20 /EHsc /W4 /utf-8 /openmp /DNDEBUG /DRAYT_SPECTRAL_MODE %INC% compile_check.cpp %GGX% /Fe:check_spectral.exe
del *.obj 2>nul
