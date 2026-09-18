"""Check the actual Fortran collision against an independent Hermite projection.

Requires numpy and an OpenACC-capable gfortran. All compiler/run artifacts use an
ASCII temporary directory. No production source or case settings are rewritten.
"""
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / '2DRBOpenacc.F90'
BASE = ROOT.parent / '均匀网格' / '2DRBOpenacc.F90'
OUT = ROOT / 'verification'
FC = shutil.which('gfortran')
if not FC:
    raise RuntimeError('gfortran is required')
ENV = dict(os.environ, ACC_DEVICE_TYPE='host')
ENV['PATH'] = str(Path(FC).parent) + os.pathsep + ENV.get('PATH', '')
WORK = Path(tempfile.mkdtemp(prefix='ptrt_verify_'))
TEXT = SOURCE.read_text(encoding='utf-8')
BASE_TEXT = BASE.read_text(encoding='utf-8')


def run(args, cwd=WORK, input=None):
    result = subprocess.run(args, cwd=cwd, input=input, text=True,
                            capture_output=True, env=ENV, timeout=120)
    if result.returncode:
        raise RuntimeError(f'{args}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def routine(text, name):
    return re.search(rf'(?im)^\s*subroutine {name}\([^\n]*\).*?'
                     rf'^\s*end subroutine {name}\b', text, re.S).group()


def compile_file(name, text, flags=()):
    src = WORK / (name + '.F90')
    src.write_text(text, encoding='utf-8')
    exe = WORK / (name + '.exe')
    run([FC, '-cpp', '-O1', '-ffree-line-length-none', '-fcheck=all',
         '-ffpe-trap=invalid,zero,overflow', *flags, str(src), '-o', str(exe)])
    return exe


E = np.array([(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),
              (1, 1), (-1, 1), (-1, -1), (1, -1)], dtype=float)
W = np.array([4/9] + [1/9]*4 + [1/36]*4)
CS = 1/3
H2 = np.einsum('ia,ib->iab', E, E) - CS*np.eye(2)
H3 = np.zeros((9, 2, 2, 2))
for a, b, c in itertools.product(range(2), repeat=3):
    H3[:, a, b, c] = E[:, a]*E[:, b]*E[:, c] - CS*(
        E[:, a]*(b == c) + E[:, b]*(a == c) + E[:, c]*(a == b))
PHI = (E[:, 0]**2-CS)*(E[:, 1]**2-CS)


def equilibrium(rho, u):
    return W*rho*(1 + E@u/CS + np.einsum('iab,a,b->i', H2, u, u)/(2*CS**2)
                  + np.einsum('iabc,a,b,c->i', H3, u, u, u)/(6*CS**3))


def projections(h):
    p1 = W*(E@(E.T@h))/CS
    a2 = np.einsum('iab,i->ab', H2, h)
    a3 = np.einsum('iabc,i->abc', H3, h)
    p2 = W*np.einsum('iab,ab->i', H2, a2)/(2*CS**2)
    p3 = W*np.einsum('iabc,abc->i', H3, a3)/(6*CS**3)
    return p1, p2, p3, a2, a3


DRIVER = '''
program verify_collision
  use commondata
  implicit none
  integer :: i,j
  allocate(f(nx,ny,0:8),f_post(0:nx+1,0:ny+1,0:8))
  allocate(rho(nx,ny),u(0:nx+1,0:ny+1),v(0:nx+1,0:ny+1),T(0:nx+1,0:ny+1))
  allocate(Fx(nx,ny),Fy(nx,ny))
  omega(0)=4.0d0/9.0d0
  omega(1:4)=1.0d0/9.0d0
  omega(5:8)=1.0d0/36.0d0
  u=0; v=0; T=0
  do j=1,ny
    do i=1,nx
      read(*,*) rho(i,j),u(i,j),v(i,j),T(i,j),f(i,j,:)
    enddo
  enddo
  !$acc enter data copyin(f,rho,u,v,T,ex,ey,opposite,omega) create(f_post,Fx,Fy)
  call collision()
  !$acc wait(1)
  !$acc update self(f_post,Fx,Fy)
  write(*,'(3es25.16)') Snu,Sq,gBeta
  do j=1,ny
    do i=1,nx
      write(*,'(11es25.16)') f_post(i,j,:),Fx(i,j),Fy(i,j)
    enddo
  enddo
  !$acc exit data delete(f,rho,u,v,T,ex,ey,opposite,omega,f_post,Fx,Fy)
end program verify_collision
'''


def collision_tests():
    # Use production commondata and collision verbatim; the oracle uses tensors,
    # not the production parity/ghost subtraction algorithm.
    module = TEXT[:TEXT.index('    end module commondata')+len('    end module commondata')]
    kernel = module + '\n' + DRIVER + '\n' + routine(TEXT, 'collision')
    rng = np.random.default_rng(260918)
    rows = []
    for case in range(128):
        rho = rng.uniform(.9, 1.1)
        u = rng.uniform(-.07, .07, 2)
        t = 0.0 if case % 4 == 0 else rng.uniform(-.5, .5)
        h = rng.normal(0, 2e-4, 9)
        h -= W*h.sum()
        if case < 2:
            h[:] = 0  # equilibrium and pure-ghost pairs
            t = 0
        f = equilibrium(rho, u) + h
        rows.append(np.r_[rho, u, t, f])
        rows.append(np.r_[rho, u, t, f + rng.uniform(-.02, .02)*W*PHI])
    data = np.array(rows)
    input_text = '\n'.join(' '.join(f'{x:.17e}' for x in r) for r in data) + '\n'
    results = []
    for ra in (1000, 1000000, 10000000000):
        for force_mode in ('buoyancy', 'magnetic'):
            # Existing optional force branch supplies nonzero Fx as well as Fy.
            extra = ['-DSideHeatedHa'] if force_mode == 'magnetic' else []
            flags = ['-DNX_OVERRIDE=16', '-DNY_OVERRIDE=16',
                     f'-DRAYLEIGH_OVERRIDE={ra}.0d0', *extra]
            outputs = []
            for backend in ('serial', 'openacc_host'):
                test_kernel = kernel.replace('(0.0d0)*(pi/180.0d0)', '(30.0d0)*(pi/180.0d0)')
                exe = compile_file(f'kernel_{ra}_{force_mode}_{backend}', test_kernel,
                                   flags + (['-fopenacc'] if backend == 'openacc_host' else []))
                text = run([str(exe)], input=input_text)
                lines = text.splitlines()
                snu, sq, gbeta = np.fromstring(lines[0], sep=' ')
                output = np.array([np.fromstring(line, sep=' ') for line in lines[1:]])
                assert output.shape == (256, 11)
                outputs.append(output)
                max_error = 0.0
                for row, got in zip(data, output):
                    rho, ux, uy, t = row[:4]
                    u, f = np.array([ux, uy]), row[4:]
                    feq = equilibrium(rho, u)
                    p1, p2, p3, a2, a3 = projections(f-feq)
                    force = got[9:]
                    viscosity = (1/snu-.5)/3
                    bmag = 400*viscosity/(16*16) if force_mode == 'magnetic' else 0
                    sine, cosine = .5, np.sqrt(3)/2
                    np.testing.assert_allclose(force, [bmag*(uy*sine*cosine-ux*sine*sine),
                        rho*gbeta*t+rho*bmag*(ux*sine*cosine-uy*cosine*cosine)], atol=2e-18)
                    a = W*(E@force)/CS
                    b = W*((E@u)*(E@force)/CS**2-u@force/CS)
                    ref = feq+p1+(1-snu)*p2+(1-sq)*p3+a+(1-snu/2)*b
                    max_error = max(max_error, float(np.max(np.abs(ref-got[:9]))))
                    np.testing.assert_allclose(got[:9], ref, rtol=0, atol=3e-15)
                    np.testing.assert_allclose(got[:9].sum(), f.sum(), rtol=0, atol=3e-15)
                    np.testing.assert_allclose(E.T@(got[:9]-f), force, rtol=0, atol=3e-15)
                    _, _, _, post2, post3 = projections(got[:9]-feq)
                    np.testing.assert_allclose(post2, (1-snu)*a2+(1-snu/2)*(
                        np.outer(u, force)+np.outer(force, u)), rtol=0, atol=3e-15)
                    np.testing.assert_allclose(post3, (1-sq)*a3, rtol=0, atol=3e-15)
                    np.testing.assert_allclose(PHI@(got[:9]-feq), 0, rtol=0, atol=3e-15)
                # Adding an arbitrary ghost mode must not change the collision output.
                np.testing.assert_allclose(output[::2], output[1::2], rtol=0, atol=3e-15)
                results.append(dict(ra=ra, force=force_mode, backend=backend,
                                    samples=256, max_tensor_error=max_error, Snu=snu, Sq=sq))
            np.testing.assert_allclose(outputs[0], outputs[1], rtol=0, atol=3e-15)
    return results


SMOKE = '''
program verify_smoke
  use commondata
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer :: k
  real(kind=8) :: mass0,mass1
  call initial()
  mass0=sum(f)
  call enter_data_2d_openacc()
  do k=1,2000
    call collision()
    call streaming()
    call bounceback()
    call macro()
    call collisionT()
    call streamingT()
    call bouncebackT()
    call macroT()
  enddo
  !$acc wait(1)
  !$acc update self(f,g,rho,u,v,T)
  if (.not.all(ieee_is_finite(f))) error stop 'nonfinite f'
  if (.not.all(ieee_is_finite(g))) error stop 'nonfinite g'
  if (.not.all(ieee_is_finite(u))) error stop 'nonfinite u'
  if (.not.all(ieee_is_finite(v))) error stop 'nonfinite v'
  if (.not.all(ieee_is_finite(T))) error stop 'nonfinite T'
  if (minval(rho).le.0.0d0) error stop 'nonpositive density'
  mass1=sum(f)
  if (abs(mass1/mass0-1.0d0).gt.1.0d-10) error stop 'mass drift'
  call calNuRe()
  write(*,'(a,7es25.16)') 'SMOKE ', mass1/mass0-1.0d0,minval(rho),maxval(rho), &
    maxval(sqrt(u*u+v*v)),minval(T),maxval(T),2000.0d0/timeUnit
  open(unit=88,file='final_fields.dat',status='replace')
  write(88,'(es25.16)') rho,u,v,T
  close(88)
  call exit_data_2d_openacc()
end program verify_smoke
'''


def smoke_tests():
    results = []
    fields = {}
    for label, source in [('baseline_luo', BASE_TEXT.replace('#define EnableUseG\n',
                        '!#define EnableUseG\n').replace('!#define EnableLegacyThermalScheme\n',
                        '#define EnableLegacyThermalScheme\n')), ('ptrt', TEXT)]:
        s = re.sub(r'(?ims)^    program main\b.*?^    end program main\b', SMOKE, source)
        for backend in ('serial', 'openacc_host'):
            flags = ['-DNX_OVERRIDE=32', '-DNY_OVERRIDE=32', '-DRAYLEIGH_OVERRIDE=10000']
            if backend == 'openacc_host':
                flags.append('-fopenacc')  # use openacc module even in serial harness
            else:
                flags.append('-lgomp')
            exe = compile_file(f'smoke_{label}_{backend}', s, flags)
            run_dir = WORK / f'run_{label}_{backend}'
            run_dir.mkdir()
            text = run([str(exe)], cwd=run_dir)
            (OUT / f'smoke_{label}_{backend}.log').write_text(text, encoding='utf-8')
            line = next(line for line in text.splitlines() if line.startswith('SMOKE '))
            values = np.fromstring(line[6:], sep=' ')
            fields[label, backend] = np.loadtxt(run_dir / 'final_fields.dat')
            results.append(dict(case=label, backend=backend, grid=[32, 32], Ra=10000,
                                steps=2000, mass_relative_drift=values[0], rho_min=values[1],
                                rho_max=values[2], speed_max=values[3], T_min=values[4],
                                T_max=values[5], t_ff=values[6]))
        np.testing.assert_allclose(fields[label, 'serial'], fields[label, 'openacc_host'],
                                   rtol=0, atol=2e-13)
    return results, float(np.max(np.abs(fields['ptrt', 'serial']-fields['baseline_luo', 'serial'])))


if __name__ == '__main__':
    # Temperature collision and thermal walls must remain byte-for-byte identical.
    for name in ('collisionT', 'streamingT', 'bouncebackT', 'macroT', 'macro', 'bounceback', 'streaming'):
        assert routine(TEXT, name) == routine(BASE_TEXT, name), name
    print('Unchanged thermal, streaming, wall, and macro routines: PASS', flush=True)
    results = dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                   baseline_sha256=hashlib.sha256(BASE.read_bytes()).hexdigest(),
                   compiler=run([FC, '--version']).splitlines()[0], work_dir=str(WORK))
    results['collision_tests'] = collision_tests()
    print('Actual Fortran collision vs Hermite tensors, conservation, ghost removal: PASS', flush=True)
    results['smoke_tests'], results['baseline_max_field_difference'] = smoke_tests()
    results['status'] = 'PASS; local serial and OpenACC host only; no GPU validation'
    (OUT / 'results.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    print(json.dumps(results, indent=2))
