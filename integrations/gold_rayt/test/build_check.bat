@echo off
rem Compiles the Material against GoLD_rayt's headers and GGX.cpp (RGB and spectral builds).
rem   GOLD_RAYT : GoLD_rayt source root (the folder that holds include\, src\ and external\). Required.
rem   VCVARS    : vcvars64.bat to use (default: Visual Studio 2022 Community).
if not defined GOLD_RAYT (
    echo set GOLD_RAYT to the GoLD_rayt source root ^(the folder with include\, src\ and external\^)
    exit /b 1
)
if not defined VCVARS set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "%VCVARS%" >nul 2>&1
cd /d %~dp0
set INC=/I ..\include /I "%GOLD_RAYT%\include" /I "%GOLD_RAYT%\external"
set GGX="%GOLD_RAYT%\src\Microfacet\GGX.cpp"
echo === RGB build
cl /nologo /O2 /std:c++20 /EHsc /W4 /utf-8 /openmp /DNDEBUG %INC% compile_check.cpp %GGX% /Fe:check_rgb.exe
echo === spectral build
cl /nologo /O2 /std:c++20 /EHsc /W4 /utf-8 /openmp /DNDEBUG /DRAYT_SPECTRAL_MODE %INC% compile_check.cpp %GGX% /Fe:check_spectral.exe
del *.obj 2>nul
