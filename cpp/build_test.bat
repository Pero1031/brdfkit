@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /O2 /EHsc /std:c++17 /W4 /utf-8 test_decoder.cpp /Fe:test_decoder.exe
