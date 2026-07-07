$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $PSScriptRoot
$sourceNames = @('2DRBOpenmpLBMCDE.F90', '2DRBOpenaccLBMCDE.F90')

$checks = @(
    @{ Name = 'temperature distribution uses D2Q9'; Pattern = 'g\(0:8' },
    @{ Name = 'temperature post-collision uses D2Q9'; Pattern = 'g_post\(0:8' },
    @{ Name = 'flow source has chi_s parameter'; Pattern = 'chi_s' },
    @{ Name = 'scalar source has chi_kappa parameter'; Pattern = 'chi_kappa' },
    @{ Name = 'flow BGK relaxation time is explicit'; Pattern = 'taufL' },
    @{ Name = 'scalar BGK relaxation time is explicit'; Pattern = 'taugL' },
    @{ Name = 'paper-style strain-rate helper exists'; Pattern = 'subroutine compute_strain_rate' },
    @{ Name = 'paper-style temperature-gradient helper exists'; Pattern = 'subroutine compute_temperature_gradient' },
    @{ Name = 'side-heated cavity case is selected'; Pattern = '#define SideHeatedCell' },
    @{ Name = 'vertical isothermal walls are selected'; Pattern = '#define VerticalWallsConstT' },
    @{ Name = 'horizontal adiabatic walls are selected'; Pattern = '#define HorizontalWallsAdiabatic' },
    @{ Name = 'side-heated temperature is centered on Tref'; Pattern = 'Thot=0\.5d0,\s*Tcold=-0\.5d0' },
    @{ Name = 'flow feq uses weakly compressible delta-rho form'; Pattern = 'deltaRhoLoc\s*=\s*rhoLoc-rho0' },
    @{ Name = 'flow macro reconstructs total rho from delta rho'; Pattern = 'rho\(i,j\)\s*=\s*rho0\s*\+\s*f\(0,i,j\)' },
    @{ Name = 'flow macro velocity divides by rho0'; Pattern = 'u\(i,j\)\s*=\s*\(momx\+0\.5d0\*Fx\(i,j\)\)/rho0' },
    @{ Name = 'Boussinesq force uses rho0 in the LBM-CDE momentum equation'; Pattern = 'Fy\(i,j\)\s*=\s*rho0\*gBeta\*\(T\(i,j\)-Tref\)' },
    @{ Name = 'isothermal thermal boundary includes lambdaT pressure term'; Pattern = '(q?lambdaT)\(1\)\*Thot\*pressure\(1,j\)/\(rho0\*cs2\)' },
    @{ Name = 'D2Q9 mixed thermal corners are explicitly repaired'; Pattern = 'D2Q9 corner diagonals belong to both' }
)

foreach ($sourceName in $sourceNames) {
    $source = Join-Path $repoRoot $sourceName
    if (-not (Test-Path -LiteralPath $source)) {
        throw "Missing $sourceName"
    }

    $text = Get-Content -LiteralPath $source -Raw

    foreach ($check in $checks) {
        if ($text -notmatch $check.Pattern) {
            throw "Static check failed in ${sourceName}: $($check.Name)"
        }
    }

    $lambdaPattern = '(?s)lambdaT\(0:8\)=\(/.*?-5\.0d0/9\.0d0.*?1\.0d0/9\.0d0.*?1\.0d0/36\.0d0.*?/\)'
    if ($text -notmatch $lambdaPattern) {
        throw "Static check failed in ${sourceName}: Chai-Zhao/LBM-CDE lambdaT weights are not present"
    }

    if ($sourceName -match 'Openacc') {
        $accChecks = @(
            @{ Name = 'OpenACC module is used'; Pattern = 'use openacc' },
            @{ Name = 'OpenACC data region is present'; Pattern = '!\$acc enter data copyin' },
            @{ Name = 'OpenACC kernels are present'; Pattern = '!\$acc parallel loop' },
            @{ Name = 'D2Q9 lookup helpers are available on device'; Pattern = 'function qomega' }
        )
        foreach ($check in $accChecks) {
            if ($text -notmatch $check.Pattern) {
                throw "Static check failed in ${sourceName}: $($check.Name)"
            }
        }
        if ($text -match '!\$omp|use omp_lib|omp_get|omp_set') {
            throw "Static check failed in ${sourceName}: OpenMP remnants found"
        }
    }
}

$ex = @(0, 1, 0, -1, 0, 1, -1, -1, 1)
$ey = @(0, 0, 1, 0, -1, 1, 1, -1, -1)
$lambdaT = @(
    (-5.0/9.0),
    (1.0/9.0), (1.0/9.0), (1.0/9.0), (1.0/9.0),
    (1.0/36.0), (1.0/36.0), (1.0/36.0), (1.0/36.0)
)

$sum0 = 0.0
$sumX = 0.0
$sumY = 0.0
$sumXX = 0.0
$sumYY = 0.0
$sumXY = 0.0
for ($i = 0; $i -lt 9; $i++) {
    $sum0 += $lambdaT[$i]
    $sumX += $lambdaT[$i] * $ex[$i]
    $sumY += $lambdaT[$i] * $ey[$i]
    $sumXX += $lambdaT[$i] * $ex[$i] * $ex[$i]
    $sumYY += $lambdaT[$i] * $ey[$i] * $ey[$i]
    $sumXY += $lambdaT[$i] * $ex[$i] * $ey[$i]
}

$tol = 1.0e-12
if ([Math]::Abs($sum0) -gt $tol -or [Math]::Abs($sumX) -gt $tol -or [Math]::Abs($sumY) -gt $tol) {
    throw "Static check failed: lambdaT zeroth/first moments are inconsistent"
}
if ([Math]::Abs($sumXX - 1.0/3.0) -gt $tol -or [Math]::Abs($sumYY - 1.0/3.0) -gt $tol -or [Math]::Abs($sumXY) -gt $tol) {
    throw "Static check failed: lambdaT second moments are inconsistent with cs2"
}

Write-Output "2DRBOpenmpLBMCDE and 2DRBOpenaccLBMCDE static checks passed"
