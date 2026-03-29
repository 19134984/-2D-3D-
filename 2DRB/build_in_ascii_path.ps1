param(
    [string]$SourceFile = "2DRB.F90",
    [string]$BuildRoot = "C:\msys64\tmp\fortran-2drb-build",
    [switch]$Run
)

$ErrorActionPreference = "Stop"

function Convert-ToMsysPath {
    param([string]$WindowsPath)

    $full = [System.IO.Path]::GetFullPath($WindowsPath)
    $drive = $full.Substring(0, 1).ToLowerInvariant()
    $rest  = $full.Substring(2).Replace("\", "/")
    return "/$drive$rest"
}

$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$sourcePath  = Join-Path $projectRoot $SourceFile

if (-not (Test-Path -LiteralPath $sourcePath)) {
    throw "Source file not found: $sourcePath"
}

if (-not (Test-Path -LiteralPath "C:\msys64\usr\bin\bash.exe")) {
    throw "MSYS2 bash was not found at C:\msys64\usr\bin\bash.exe"
}

if (-not (Test-Path -LiteralPath "C:\msys64\ucrt64\bin\gfortran.exe")) {
    throw "gfortran was not found at C:\msys64\ucrt64\bin\gfortran.exe"
}

New-Item -ItemType Directory -Force -Path $BuildRoot | Out-Null

$sourceName = Split-Path -Leaf $sourcePath
$exeName    = ([System.IO.Path]::GetFileNameWithoutExtension($sourceName)) + ".exe"
$buildSrc   = Join-Path $BuildRoot $sourceName

Copy-Item -LiteralPath $sourcePath -Destination $buildSrc -Force

$buildRootMsys = Convert-ToMsysPath $BuildRoot
$compileCmd = "export PATH=/ucrt64/bin:/usr/bin:`$PATH && cd '$buildRootMsys' && /ucrt64/bin/gfortran -cpp -fopenmp '$sourceName' -o '$exeName'"

if ($Run) {
    $compileCmd += " && './$exeName'"
}

Write-Host "Build root  : $BuildRoot"
Write-Host "Source file : $sourcePath"
Write-Host "Executable  : $(Join-Path $BuildRoot $exeName)"

& "C:\msys64\usr\bin\bash.exe" -lc $compileCmd
