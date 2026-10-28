"""Actual-source OpenACC host checks for both closed-cavity conduction geometries.
Adapted from the independent review driver; overrides buoyancy and the initial RB perturbation.
Outputs stay under _local_checks/conduction; never modifies the production source.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess

from verify_transport import replace_parameter, COMPILER, ENV

ROOT = Path(__file__).resolve().parent
HERE = ROOT / '_local_checks' / 'conduction'
SOURCE = ROOT / '2DRBOpenaccEntropic.F90'

DRIVER = r'''
program main
  use openacc
  use commondata
  implicit none
  real(8) :: initialTemperature(nx,ny),initialMass,maxTemperatureError,maxVelocity,maxMassError
  integer :: n
  call acc_init(acc_device_host)
  call initial()
  initialTemperature=T
  initialMass=sum(rho)
  call enter_data_2d_openacc()
  maxTemperatureError=0d0; maxVelocity=0d0; maxMassError=0d0
  do n=1,600
    call collision()
    call collisionT()
    call streaming()
    call bounceback()
    call streamingT()
    call bouncebackT()
    call macro()
    call update_host_snapshot_2d_openacc()
    maxTemperatureError=max(maxTemperatureError,maxval(abs(T-initialTemperature)))
    maxVelocity=max(maxVelocity,maxval(abs(u)),maxval(abs(v)))
    maxMassError=max(maxMassError,abs(sum(rho)-initialMass))
  enddo
  call check_fields_entropic()
#ifdef SideHeatedCell
  call SideHeatedcalc_Nu_global()
  call SideHeatedcalc_Nu_wall_avg()
#else
  call RBcalc_Nu_global()
  call RBcalc_Nu_wall_avg()
#endif
  write(*,'(A,7ES25.16E3)') 'CONDUCTION_RESULT ',maxTemperatureError,maxVelocity,maxMassError, &
      Nu_global,Nu_hot,Nu_cold,Nu_middle
  if(maxTemperatureError>1d-11.or.maxVelocity>1d-11.or.maxMassError>1d-9) &
      error stop 'Conduction profile, velocity, or mass regression'
  if(maxval(abs((/Nu_global,Nu_hot,Nu_cold,Nu_middle/)-1d0))>1d-9) &
      error stop 'Conduction Nu regression'
  call exit_data_2d_openacc()
end program main
'''


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    raw = SOURCE.read_bytes()
    text = raw.decode('utf-8-sig').replace('\r\n', '\n')
    text, count = re.subn(r'(?ims)^\s*program main\s*$.*?^\s*end program main\s*$', DRIVER, text)
    assert count == 1
    text = re.sub(r'(?im)^\s*integer\(kind=4\),\s*parameter\s*::\s*nx=\d+,\s*ny=\d+.*$',
                  '        integer(kind=4), parameter :: nx=16, ny=12', text)
    for name, value in {'gBeta': '0d0', 'timeUnit': '1d0', 'velocityUnit': '1d0',
                        'rbInitialPerturbation': '0d0', 'tauf': '0.8d0'}.items():
        text = replace_parameter(text, name, value)
    cases = []
    for geometry in ['SideHeatedCell', 'RayleighBenardCell']:
        case_text = text
        if geometry == 'RayleighBenardCell':
            case_text = case_text.replace('#define SideHeatedCell', '!#define SideHeatedCell', 1)
            case_text = case_text.replace('!#define RayleighBenardCell', '#define RayleighBenardCell', 1)
        folder = HERE / ('entropic_conduction_' + geometry)
        folder.mkdir(exist_ok=True)
        source = folder / 'conduction.F90'
        source.write_text(case_text, encoding='utf-8')
        compile_result = subprocess.run([str(COMPILER), '-cpp', '-ffree-line-length-none', '-fopenacc',
                                         '-fcheck=all', '-fbacktrace', '-O2', source.name, '-o', 'conduction.exe'],
                                        cwd=folder, env=ENV, capture_output=True, text=True)
        if compile_result.returncode:
            raise RuntimeError(compile_result.stdout + compile_result.stderr)
        env = dict(ENV, ACC_DEVICE_TYPE='host')
        run = subprocess.run([str(folder / 'conduction.exe')], cwd=folder, env=env,
                             capture_output=True, text=True, timeout=60)
        (folder / 'console.txt').write_text(run.stdout + run.stderr, encoding='utf-8')
        if run.returncode:
            raise RuntimeError(run.stdout + run.stderr)
        values = next(line.split()[1:] for line in run.stdout.splitlines() if line.startswith('CONDUCTION_RESULT'))
        names = ['maximum_temperature_error', 'maximum_absolute_velocity', 'maximum_absolute_mass_error',
                 'Nu_global', 'Nu_hot', 'Nu_cold', 'Nu_middle']
        cases.append(dict(geometry=geometry, nx=16, ny=12, steps=600, backend='GNU OpenACC host',
                          **dict(zip(names, map(float, values)))))
    result = dict(status='passed', source_sha256=hashlib.sha256(raw).hexdigest(), cases=cases)
    (ROOT / 'conduction_validation_local.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
