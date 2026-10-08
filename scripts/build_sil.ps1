param([string]$VsRoot='C:/Program Files (x86)/Microsoft Visual Studio/2019/BuildTools')
$ErrorActionPreference='Stop'
$ProjectRoot=Split-Path $PSScriptRoot -Parent
$SilBuild=Join-Path $ProjectRoot 'validation/sil'
New-Item -ItemType Directory -Path $SilBuild -Force | Out-Null
$CompilerInit=Join-Path $VsRoot 'VC/Auxiliary/Build/vcvars64.bat'
if (!(Test-Path -LiteralPath $CompilerInit)) { throw 'MSVC x64 environment not found; pass -VsRoot.' }
$Source=Join-Path $PSScriptRoot 'sil_core.c'
$PlantSource=Join-Path $PSScriptRoot 'sil_plant.c'
$Dll=Join-Path $SilBuild 'sil_core.dll'
$SilObjects=$SilBuild.Replace('\','/')+'/'
$Batch=Join-Path $SilBuild 'compile.cmd'
$Body="@echo off`r`ncall `"$CompilerInit`"`r`nif errorlevel 1 exit /b 1`r`ncl /nologo /O2 /fp:precise /W4 /LD `"$Source`" `"$PlantSource`" /Fo:`"$SilObjects`" /link /OUT:`"$Dll`"`r`nexit /b %errorlevel%`r`n"
Set-Content -LiteralPath $Batch -Value $Body -Encoding Ascii
& $env:ComSpec /d /c $Batch 2>&1 | Tee-Object -FilePath (Join-Path $SilBuild 'compile.log')
if($LASTEXITCODE -ne 0){throw 'SIL compile failed.'}
