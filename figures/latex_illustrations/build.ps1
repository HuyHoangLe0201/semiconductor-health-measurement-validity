param([string]$TexBin = 'D:/texlive/2026/bin/windows')
$ErrorActionPreference = 'Stop'
$latex = Join-Path $TexBin 'pdflatex.exe'
if (-not (Test-Path -LiteralPath $latex)) {
    $latex = (Get-Command pdflatex -ErrorAction Stop).Source
}
Push-Location $PSScriptRoot
try {
    1..2 | ForEach-Object {
        & $latex '-interaction=nonstopmode' '-halt-on-error' 'measurement_illustrations.tex'
        if ($LASTEXITCODE -ne 0) { throw 'LaTeX compilation failed; inspect measurement_illustrations.log.' }
    }
    Write-Output (Join-Path $PSScriptRoot 'measurement_illustrations.pdf')
} finally { Pop-Location }
