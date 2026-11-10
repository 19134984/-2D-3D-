"""验证新HRR副本：已有公式、双级退化对照、三级导热、输出间隔和重启。"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import re
import shutil
import struct
import subprocess
import sys

import numpy as np

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'运行脚本/run_cases/Multilevel_HRR0995_Ra9_Ra10_N2048_20261103'
WORK=Path('C:/Users/User/.codex/visualizations/2026/10/26/01a176f3-942f-7de2-84f7-0618b0d3d607/multilevel_hrr_20261103')
WORK.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('existing_multilevel_verifier',ROOT/'运行脚本/tools/verify_multilevel.py')
v=importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
manifest=json.loads((OUT/'manifest.json').read_text(encoding='utf-8'))
v.SOURCE=(ROOT/manifest['new_solver']).read_text(encoding='utf-8')
v.PARENT=(ROOT/manifest['hrr_parent']).read_text(encoding='utf-8')
FC='C:/msys64/ucrt64/bin/gfortran.exe'
ENV=dict(os.environ,ACC_DEVICE_TYPE='host',OMP_NUM_THREADS='1')
ENV['PATH']=str(Path(FC).parent)+os.pathsep+ENV['PATH']
REPORT=dict(backend='GCC OpenACC host; NVIDIA device checks are separate',build_root=str(WORK),tests={})


def save():
    (OUT/'local_validation.json').write_text(json.dumps(REPORT,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def run(args,folder,expect=0,log='run.log'):
    p=subprocess.run(args,cwd=folder,env=ENV,capture_output=True,timeout=180)
    (folder/log).write_bytes(p.stdout+p.stderr)
    if expect==0: assert p.returncode==0,(folder,(p.stdout+p.stderr).decode(errors='replace')[-3500:])
    else: assert p.returncode!=0,(folder,'Expected rejection')
    return p


def build(folder,text,program=None):
    folder.mkdir(parents=True,exist_ok=True)
    (folder/'solver.F90').write_text(text,encoding='utf-8')
    flags=['-cpp','-DREG_HRR','-fopenacc','-ffree-line-length-none','-O0','-g',
           '-fcheck=all','-ffpe-trap=invalid,zero,overflow']
    if program:
        library,count=re.subn(r'^    program main\b.*?^    end program main\b','',text,flags=re.M|re.S)
        assert count==1
        (folder/'library.F90').write_text(library,encoding='utf-8')
        (folder/'probe.F90').write_text(program,encoding='utf-8')
        sources=['library.F90','probe.F90']
    else: sources=['solver.F90']
    run([FC,*flags,'-fsyntax-only',*sources],folder,log='syntax.log')
    run([FC,*flags,*sources,'-o','solver.exe'],folder,log='compile.log')


def variant(levels=3,n=128,fine=16,inner=8,steps=32,old=False,**kwargs):
    text=v.variant(levels=levels,n=n,fine=fine,inner=inner,steps=steps,old=old,**kwargs)
    text=re.sub(r'\bNuReSampleInterval=[\d.]+d0','NuReSampleInterval=8.0d0/timeUnit',text)
    text=re.sub(r'\bRayleigh=[\d.]+d[\d]+','Rayleigh=1.0d6',text,count=1)
    return text


def execute(label,text,expect=0):
    folder=WORK/label
    build(folder,text)
    p=run([str(folder/'solver.exe')],folder,expect=expect)
    if expect==0:
        final=folder/'FinalNuRe_2DOpenaccMultiblock.dat'
        values=np.loadtxt(final) if final.exists() else np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat')[-1,2:]
        assert np.isfinite(values).all()
        REPORT['tests'][label]=dict(returncode=0,final_metrics=values.tolist())
    else:
        REPORT['tests'][label]=dict(returncode=p.returncode,rejection=p.stdout.decode(errors='replace')[-500:])
    save()
    print(label,'PASS',flush=True)
    return folder


def checkpoint(folder):
    return folder/(folder/'reloadFile2DOpenaccMultilevelHRR-latest.meta').read_text().splitlines()[1]


def inspect_checkpoint(path):
    raw=path.read_bytes()
    assert raw[:16].rstrip()==b'MLHRRESTART0001'
    config=struct.unpack_from('<16i',raw,16)
    offset=16+64+96
    total,samples=struct.unpack_from('<2i',raw,offset); offset+=8
    offset+=16
    nb,sigma=struct.unpack_from('<id',raw,offset); offset+=12
    assert nb==4*(config[2]-1)+1 and sigma==.995
    shapes=[]
    for _ in range(nb):
        level,ni,nj,gi,gj=struct.unpack_from('<5i',raw,offset); offset+=20
        dx,x,y=struct.unpack_from('<3d',raw,offset); offset+=24
        count=ni*nj*20
        data=np.frombuffer(raw,dtype='<f8',count=count,offset=offset)
        assert np.isfinite(data).all()
        offset+=8*count
        shapes.append([level,ni,nj,dx])
    for l in range(1,config[2]):
        nf,cycles=struct.unpack_from('<2i',raw,offset); offset+=8
        assert cycles==total//2**l
        count=nf*20*3
        assert np.isfinite(np.frombuffer(raw,dtype='<f8',count=count,offset=offset)).all()
        offset+=8*count
    assert offset==len(raw)
    return dict(total_steps=total,samples=samples,sigma=sigma,blocks=shapes,all_fields_and_interface_history_finite=True)


def main():
    # 与已有独立Hermite参考数据逐节点比较，包括非零当前力、速度梯度与接口回环。
    folder=WORK/'kernel_formula'
    for name in ['input.dat','expected.dat']:
        folder.mkdir(exist_ok=True)
        shutil.copyfile(OUT/'verification'/name,folder/name)
    build(folder,v.SOURCE,(OUT/'verification/verify.F90').read_text(encoding='utf-8'))
    run([str(folder/'solver.exe')],folder)
    error=float((folder/'verification_error.dat').read_text())
    assert error<1e-12
    REPORT['tests']['kernel_formula']=dict(node_count=252,max_absolute_error=error)
    print('kernel_formula PASS',error,flush=True); save()
    for c in manifest['cases']:
        text=(OUT/c['case']/'source/solver.F90').read_text(encoding='utf-8')
        folder=WORK/c['case']
        probe='''program parameters
use multilevel_solver
implicit none
integer :: b,l
call initial_geometry()
write(*,*) nx,ny,nGridLevels,nBlocks,loadInitField,itc_max,NuReSampleIntervalItc,outputSnapshotIntervalItc
write(*,*) Rayleigh,Prandtl,Mach,hrrMixingSigma,timeUnit,viscosity,diffusivity,paraA
do l=1,nGridLevels
write(*,*) grid(l)%dx,SnuLevel(l),SqLevel(l),QkLevel(l),QnuLevel(l)
enddo
do b=1,nBlocks
write(*,*) block(b)%level,block(b)%ni,block(b)%nj
enddo
end program
'''
        build(folder,text,probe)
        result=run([str(folder/'solver.exe')],folder)
        lines=result.stdout.decode().splitlines()
        ints=[int(x) for x in lines[0].split()]
        physics=[float(x) for x in lines[1].split()]
        levels=[[float(x) for x in l.split()] for l in lines[2:5]]
        shapes=[[int(x) for x in l.split()] for l in lines[5:14]]
        assert ints[:5]==[2048,2048,3,9,0]
        assert physics[:4]==[float(c['Ra']),.7,.1,.995]
        assert np.allclose([p[0] for p in levels],[1,2,4],rtol=0,atol=0)
        assert all(p[1]==p[2] for p in levels)
        assert abs(ints[6]/physics[4]-.5)<4/physics[4]
        assert abs(ints[7]/physics[4]-5)<4/physics[4]
        assert 1000<=ints[5]/physics[4]<1000.6
        stored=sum(ni*nj for _,ni,nj in shapes)
        c.update(actual_target_tff=ints[5]/physics[4],target_fine_steps=ints[5],
                 timeUnit=physics[4],block_shapes=shapes,stored_nodes=stored,
                 level_parameters=levels,
                 estimated_snapshot_payload_bytes=stored*32*201,
                 estimated_restart_field_payload_bytes=stored*160*11)
        REPORT['tests'][c['case']]=dict(syntax_and_geometry_exit=0,parameters=physics,level_parameters=levels,
            mesh_and_runtime=ints,stored_nodes=stored)
        print(c['case'],'PARAMETERS PASS',flush=True); save()
    (OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    # 退化为一/两级时，与原HRR全状态和已有Nu/Re定义比较。
    for levels,n,fine in [(1,32,8),(2,64,8)]:
        ref=execute(f'twolevel_parent_L{levels}',variant(levels=levels,n=n,fine=fine,old=True))
        new=execute(f'new_L{levels}',variant(levels=levels,n=n,fine=fine))
        a=np.fromfile(ref/'raw-state.bin',dtype='<f8')
        b=np.fromfile(new/'raw-state.bin',dtype='<f8')
        np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-13)
        ha=np.loadtxt(ref/'NuRe_2DOpenaccMultiblock.dat')
        hb=np.loadtxt(new/'NuRe_2DOpenaccMultiblock.dat')
        np.testing.assert_allclose(ha,hb,rtol=1e-10,atol=1e-11)
        REPORT['tests'][f'HRR_parent_comparison_L{levels}']=dict(full_state_max_abs_difference=float(np.max(np.abs(a-b))),
            history_max_abs_difference=float(np.max(np.abs(ha-hb))))
        print('parent comparison',levels,'PASS',flush=True); save()
    full=execute('levels3',variant())
    REPORT['tests']['levels3']['checkpoint']=inspect_checkpoint(checkpoint(full)); save()
    conduction=execute('conduction3',variant(conduction=True))
    values=np.loadtxt(conduction/'FinalNuRe_2DOpenaccMultiblock.dat')
    assert max(abs(values[[1,3,4,5]]-1))<1e-8
    REPORT['tests']['conduction3']['NuVol_hot_cold_middle']=values[[1,3,4,5]].tolist(); save()
    restart=execute('restart3',variant(steps=16))
    build(restart,variant(restart=True))
    run([str(restart/'solver.exe')],restart)
    assert checkpoint(restart).read_bytes()==checkpoint(full).read_bytes()
    np.testing.assert_array_equal(np.loadtxt(restart/'NuRe_2DOpenaccMultiblock.dat'),
                                  np.loadtxt(full/'NuRe_2DOpenaccMultiblock.dat'))
    REPORT['tests']['restart3']['continuation']='Bitwise identical checkpoint and equal history'; save()
    text=variant(steps=64)
    text=re.sub(r'\boutputSnapshotInterval=8.0d0/timeUnit', 'outputSnapshotInterval=24.0d0/timeUnit',text)
    separated=execute('separated_sampling',text)
    h=np.loadtxt(separated/'NuRe_2DOpenaccMultiblock.dat')
    np.testing.assert_array_equal(h[:,1],np.arange(8,65,8))
    snapshots=list(separated.glob('buoyancyCavity2DOpenaccMultiblockSnapshot*.bin'))
    assert len(snapshots)==3
    REPORT['tests']['separated_sampling'].update(history_steps=h[:,1].tolist(),snapshot_count=3); save()
    # 不同sigma的检查点必须拒绝，旧未正则化检查点也由独立magic拒绝。
    bad=WORK/'restart_wrong_sigma'; bad.mkdir(exist_ok=True)
    for name in [checkpoint(full).name,'reloadFile2DOpenaccMultilevelHRR-latest.meta',
                 'NuRe_2DOpenaccMultiblock.dat','OutputFiles_2DOpenaccMultiblock.dat']:
        shutil.copyfile(full/name,bad/name)
    text=variant(steps=40,restart=True).replace('hrrMixingSigma=0.995d0','hrrMixingSigma=0.99d0')
    build(bad,text); rejected=run([str(bad/'solver.exe')],bad,expect=1)
    assert b'Checkpoint HRR sigma mismatch' in rejected.stderr+rejected.stdout
    REPORT['tests']['restart_wrong_sigma']=dict(rejected=True); save()
    original=json.loads((OUT/'protected_originals_before.json').read_text(encoding='utf-8'))
    actual={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in original}
    assert actual==original
    REPORT.update(status='passed',protected_originals_unchanged=True,
        scope='Implementation regression only; long-time Ra9/Ra10 and benchmark agreement require server production runs')
    save(); print('ALL LOCAL CHECKS PASSED',flush=True)


if __name__=='__main__':
    try: main()
    except Exception as e:
        REPORT.update(status='failed',failure=str(e)); save(); raise
