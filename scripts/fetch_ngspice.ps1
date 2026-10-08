# Portable official ngspice 47 runtime; no registry or global PATH changes.
$ErrorActionPreference = 'Stop'
$PaperRoot = Split-Path $PSScriptRoot -Parent
$ExternalRoot = Join-Path $PaperRoot 'external'
$ArchivePath = Join-Path $ExternalRoot 'ngspice-47_64.7z'
$ExtractRoot = Join-Path $ExternalRoot 'ngspice47'
$ExpectedHash = '59225971bd68cdd1199443649aa4615a9e6d684933f205ab49006a3942518f5a'
New-Item -ItemType Directory -Path $ExternalRoot -Force | Out-Null
if (-not (Test-Path -LiteralPath $ArchivePath)) {
    Invoke-WebRequest -Uri 'https://sourceforge.net/projects/ngspice/files/ng-spice-rework/47/ngspice-47_64.7z/download' -OutFile $ArchivePath
}
$ActualHash = (Get-FileHash -LiteralPath $ArchivePath -Algorithm SHA256).Hash.ToLowerInvariant()
if ($ActualHash -ne $ExpectedHash) { throw 'ngspice archive hash does not match the pinned official release.' }
New-Item -ItemType Directory -Path $ExtractRoot -Force | Out-Null
tar -xf $ArchivePath -C $ExtractRoot
if ($LASTEXITCODE -ne 0) { throw 'ngspice archive extraction failed.' }
& (Join-Path $ExtractRoot 'Spice64/bin/ngspice_con.exe') --version
