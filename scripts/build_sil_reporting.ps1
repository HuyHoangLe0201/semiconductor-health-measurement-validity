param([string]$VsRoot='C:/Program Files (x86)/Microsoft Visual Studio/2019/BuildTools')
$ErrorActionPreference='Stop'
$ReportBuild=Join-Path (Split-Path $PSScriptRoot -Parent) 'validation/sil_reporting'
New-Item -ItemType Directory -Path $ReportBuild -Force | Out-Null
$ReportInit=Join-Path $VsRoot 'VC/Auxiliary/Build/vcvars64.bat'
if (!(Test-Path -LiteralPath $ReportInit)) { throw 'MSVC x64 environment not found.' }
$ReportSources=@((Join-Path $PSScriptRoot 'sil_mc_native.c'),(Join-Path $PSScriptRoot 'sil_plant.c'))
$ReportDll=Join-Path $ReportBuild 'sil_reporting.dll'
$ReportObjects=$ReportBuild.Replace('\','/')+'/'
$ReportBatch=Join-Path $ReportBuild 'compile.cmd'
$ReportBody="@echo off`r`ncall `"$ReportInit`"`r`nif errorlevel 1 exit /b 1`r`ncl /nologo /O2 /fp:precise /W4 /LD `"$($ReportSources[0])`" `"$($ReportSources[1])`" /Fo:`"$ReportObjects`" /link /OUT:`"$ReportDll`"`r`nexit /b %errorlevel%`r`n"
Set-Content -LiteralPath $ReportBatch -Value $ReportBody -Encoding Ascii
& $env:ComSpec /d /c $ReportBatch 2>&1 | Tee-Object -FilePath (Join-Path $ReportBuild 'compile.log')
if ($LASTEXITCODE -ne 0) { throw 'Reporting SIL compile failed.' }
