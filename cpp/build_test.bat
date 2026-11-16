@echo off
rem Builds the cross-check driver with MSVC. Set VCVARS to your vcvars64.bat if Visual Studio is not the 2022 Community edition.
if not defined VCVARS set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "%VCVARS%" >nul
cd /d %~dp0
cl /nologo /O2 /EHsc /std:c++17 /W4 /utf-8 test_decoder.cpp /Fe:test_decoder.exe
