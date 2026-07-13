$ErrorActionPreference = 'Stop'

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..\..')
$outputDir = Join-Path $repoRoot 'output\pdf'
$texSource = 'docs/lbm-cde-trt-derivation.tex'
$assembler = Join-Path $PSScriptRoot 'assemble_report.py'
$python = (Get-Command python -ErrorAction Stop).Source
$tempBase = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
$tempRoot = [IO.Path]::GetFullPath(
    (Join-Path $tempBase 'lbm-cde-trt-report-build')
)

if (-not $tempRoot.StartsWith($tempBase, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to use a build directory outside the system temp root: $tempRoot"
}

& $python $assembler
if ($LASTEXITCODE -ne 0) {
    throw "Markdown report assembly failed with exit code $LASTEXITCODE"
}

$env:SOURCE_DATE_EPOCH = '1783900800'
$env:FORCE_SOURCE_DATE = '1'

if (Test-Path -LiteralPath $tempRoot) {
    Remove-Item -LiteralPath $tempRoot -Recurse -Force
}
New-Item -ItemType Directory -Force (Join-Path $tempRoot 'output\pdf') | Out-Null
Copy-Item -LiteralPath (Join-Path $repoRoot 'docs') `
    -Destination (Join-Path $tempRoot 'docs') -Recurse -Force

New-Item -ItemType Directory -Force $outputDir | Out-Null

Push-Location $tempRoot
try {
    foreach ($pass in 1..2) {
        & xelatex `
            -shell-escape `
            -interaction=nonstopmode `
            -halt-on-error `
            -file-line-error `
            -jobname=lbm-cde-trt-derivation `
            -output-directory=output/pdf `
            $texSource
        if ($LASTEXITCODE -ne 0) {
            throw "XeLaTeX pass $pass failed with exit code $LASTEXITCODE"
        }
    }
}
finally {
    Pop-Location
}

$pdf = Join-Path $outputDir 'lbm-cde-trt-derivation.pdf'
$builtPdf = Join-Path $tempRoot 'output\pdf\lbm-cde-trt-derivation.pdf'
if (-not (Test-Path -LiteralPath $builtPdf)) {
    throw "Expected PDF was not created: $builtPdf"
}
Copy-Item -LiteralPath $builtPdf -Destination $pdf -Force

Get-Item -LiteralPath $pdf | Select-Object FullName, Length, LastWriteTime
