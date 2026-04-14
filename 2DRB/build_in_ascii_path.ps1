param(
    [string]$SourceFile = "2DRB.F90",
    [string]$BuildRoot = "C:\msys64\tmp\fortran-2drb-build",
    [switch]$Run,
    [switch]$SyntaxOnly,
    [switch]$UseMpi,
    [string]$ParallelFlag = "-fopenmp",
    [string[]]$ExtraArgs = @()
)

$ErrorActionPreference = "Stop"

$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$sourcePath  = Join-Path $projectRoot $SourceFile

if (-not (Test-Path -LiteralPath $sourcePath)) {
    throw "Source file not found: $sourcePath"
}

$msysBin = "C:\msys64\ucrt64\bin"
$compilerName = if ($UseMpi) { "mpifort.exe" } else { "gfortran.exe" }
$compilerPath = Join-Path $msysBin $compilerName

if (-not (Test-Path -LiteralPath $compilerPath)) {
    throw "$compilerName was not found at $compilerPath"
}

New-Item -ItemType Directory -Force -Path $BuildRoot | Out-Null

$sourceName = Split-Path -Leaf $sourcePath
$exeName    = ([System.IO.Path]::GetFileNameWithoutExtension($sourceName)) + ".exe"
$buildSrc   = Join-Path $BuildRoot $sourceName

Copy-Item -LiteralPath $sourcePath -Destination $buildSrc -Force

$env:PATH = "$msysBin;$env:PATH"
$compileArgs = @("-cpp")
if ($ParallelFlag) {
    $compileArgs += $ParallelFlag
}
if ($ExtraArgs.Count -gt 0) {
    $compileArgs += $ExtraArgs
}
if ($SyntaxOnly) {
    $compileArgs += @("-fsyntax-only", $buildSrc)
} else {
    $compileArgs += @($buildSrc, "-o", (Join-Path $BuildRoot $exeName))
}

Write-Host "Build root  : $BuildRoot"
Write-Host "Source file : $sourcePath"
Write-Host "Compiler    : $compilerPath"
if ($SyntaxOnly) {
    Write-Host "Mode        : syntax-only"
} else {
    Write-Host "Executable  : $(Join-Path $BuildRoot $exeName)"
}

& $compilerPath @compileArgs
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

if ($Run) {
    if ($SyntaxOnly) {
        throw "-Run cannot be used together with -SyntaxOnly"
    }
    & (Join-Path $BuildRoot $exeName)
    exit $LASTEXITCODE
}
