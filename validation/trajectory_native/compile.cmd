@echo off
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul
cl /nologo /O2 /fp:precise /W4 /LD tv_driver.c tv_plant.c /link /OUT:trajectory.dll
