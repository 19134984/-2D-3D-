param(
    [string]$D2Q9SourceFile = "Xs\D2Q9TRT_D2Q9BGK\2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9BGK.F90",
    [string]$D2Q5SourceFile = "Xs\D2Q9TRT_D2Q5LuoTRT\2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90"
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$compiler = "C:\msys64\ucrt64\bin\gfortran.exe"

foreach ($sourceFile in @($D2Q9SourceFile, $D2Q5SourceFile)) {
    $sourcePath = Join-Path $repoRoot $sourceFile
    if (-not (Test-Path -LiteralPath $sourcePath)) {
        throw "Source file not found: $sourcePath"
    }
}
if (-not (Test-Path -LiteralPath $compiler)) {
    throw "gfortran not found: $compiler"
}

$systemTemp = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
$testRoot = [IO.Path]::GetFullPath(
    (Join-Path $systemTemp ("lbmcde-lattice-allocation-" + [guid]::NewGuid().ToString("N")))
)
if (-not $testRoot.StartsWith($systemTemp, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to use a test directory outside the system temp root: $testRoot"
}

function Invoke-LatticeCase {
    param(
        [string]$Name,
        [string]$SourceFile,
        [string]$FlowMacro,
        [string]$ExpectedLattice,
        [string]$ExpectedFlowPolicy
    )

    $sourcePath = Join-Path $repoRoot $SourceFile
    $caseDir = Join-Path $testRoot $Name
    New-Item -ItemType Directory -Force -Path $caseDir | Out-Null
    $caseSource = Join-Path $caseDir "2DRBOpenaccLBMCDE.F90"
    $caseExe = Join-Path $caseDir "solver.exe"
    Copy-Item -LiteralPath $sourcePath -Destination $caseSource -Force

    $compileArgs = @(
        "-cpp", "-O0", "-g", "-fopenacc", "-fcheck=all", "-fbacktrace",
        "-DsteadyFlow", "-DNX_OVERRIDE=32", "-DNY_OVERRIDE=32", "-DITC_MAX_OVERRIDE=4",
        "-DCHECK_INTERVAL_OVERRIDE=2", "-DRAYLEIGH_OVERRIDE=1.0d4", $FlowMacro,
        $caseSource, "-o", $caseExe
    )

    Push-Location $caseDir
    try {
        & $compiler @compileArgs
        if ($LASTEXITCODE -ne 0) {
            throw "$Name compile failed with exit code $LASTEXITCODE"
        }

        $oldDeviceType = $env:ACC_DEVICE_TYPE
        try {
            $env:ACC_DEVICE_TYPE = "host"
            & $caseExe
            if ($LASTEXITCODE -ne 0) {
                throw "$Name bounds-check smoke run failed with exit code $LASTEXITCODE"
            }
        } finally {
            if ($null -eq $oldDeviceType) {
                Remove-Item Env:ACC_DEVICE_TYPE -ErrorAction SilentlyContinue
            } else {
                $env:ACC_DEVICE_TYPE = $oldDeviceType
            }
        }

        $settings = Get-Content -Raw -LiteralPath (Join-Path $caseDir "SimulationSettings2DOpenaccLBMCDE.txt")
        if ($settings -notmatch [regex]::Escape($ExpectedLattice)) {
            throw "$Name did not report the expected lattice: $ExpectedLattice"
        }
        if ($settings -notmatch [regex]::Escape($ExpectedFlowPolicy)) {
            throw "$Name did not report the expected flow magic policy: $ExpectedFlowPolicy"
        }
        if ($settings -notmatch "Final itc\s*=\s*4") {
            throw "$Name did not complete the four-step smoke run"
        }
    } finally {
        Pop-Location
    }
}

New-Item -ItemType Directory -Force -Path $testRoot | Out-Null
try {
    $basePolicy = "Flow odd policy = magic parameter on base/raw even collision time"
    $effectivePolicy = "Flow odd policy = magic parameter on chi_s-corrected effective even scale"
    Invoke-LatticeCase -Name "d2q9_base" -SourceFile $D2Q9SourceFile -FlowMacro "-DFLOW_ODD_BASE_MAGIC" `
        -ExpectedLattice "Temperature lattice = D2Q9 BGK LBM-CDE" -ExpectedFlowPolicy $basePolicy
    Invoke-LatticeCase -Name "d2q9_effective" -SourceFile $D2Q9SourceFile -FlowMacro "-DFLOW_ODD_EFFECTIVE_MAGIC" `
        -ExpectedLattice "Temperature lattice = D2Q9 BGK LBM-CDE" -ExpectedFlowPolicy $effectivePolicy
    Invoke-LatticeCase -Name "d2q5_base" -SourceFile $D2Q5SourceFile -FlowMacro "-DFLOW_ODD_BASE_MAGIC" `
        -ExpectedLattice "Temperature lattice = D2Q5 TRT Luo/Wang" -ExpectedFlowPolicy $basePolicy
    Invoke-LatticeCase -Name "d2q5_effective" -SourceFile $D2Q5SourceFile -FlowMacro "-DFLOW_ODD_EFFECTIVE_MAGIC" `
        -ExpectedLattice "Temperature lattice = D2Q5 TRT Luo/Wang" -ExpectedFlowPolicy $effectivePolicy
    Write-Output "PASS: both temperature solvers and both flow-magic policies completed bounds-checked OpenACC host smoke runs."
} finally {
    if (Test-Path -LiteralPath $testRoot) {
        $resolved = [IO.Path]::GetFullPath($testRoot)
        if (-not $resolved.StartsWith($systemTemp, [StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to remove a path outside the system temp root: $resolved"
        }
        Remove-Item -LiteralPath $resolved -Recurse -Force
    }
}
