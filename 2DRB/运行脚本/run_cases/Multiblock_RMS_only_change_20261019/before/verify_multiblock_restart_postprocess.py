"""Multiblock file/restart/periodic regression on small CPU cases.
Requires Python + numpy and gfortran. Uses a fresh temporary directory, never real run data.
Use --before <pre-edit F90> for a numerical comparison to the pre-edit solver.
OpenACC directives are syntax checked; these tests do not establish NVIDIA GPU correctness.
"""
from pathlib import Path
import argparse, json, os, re, shutil, struct, subprocess, tempfile
import numpy as np
parser=argparse.ArgumentParser()
parser.add_argument('--before',type=Path)
parser.add_argument('--compiler',default=shutil.which('gfortran') or r'C:\msys64\ucrt64\bin\gfortran.exe')
args=parser.parse_args()
repo=Path(__file__).resolve().parents[2]
os.chdir(repo)
root=Path(tempfile.mkdtemp(prefix='mb-regression-'))
compiler=args.compiler
fc=compiler
env=os.environ.copy()
env['PATH']=str(Path(compiler).parent)+os.pathsep+env['PATH']
src=(repo/'多块网格/2DRBOpenaccMultiblock.F90').read_text(encoding='utf-8')
print('Validation directory:',root,flush=True)
def small(s,steady=False,periodic=False,ratio=2,target=60):
 s=re.sub(r'\bnx=1024', 'nx=32',s);s=re.sub(r'\bny=1024','ny=32',s)
 for key,val in [('fineLayerCellsLeft',8),('fineLayerCellsRight',9),('fineLayerCellsBottom',8),('fineLayerCellsTop',9),('refineRatio',ratio)]:
  s=re.sub(r'\b'+key+r'=\d+',key+'='+str(val),s,count=1)
 for key,val in [('outputSnapshotInterval',6),('reloadFileInterval',20),('outputPltFileInterval',12),('unsteadyRunDuration',target)]:
  s=re.sub(r'(?m)^(\s*real\(kind=8\), parameter :: '+key+r'=)[^!\n]*',r'\g<1>'+str(val)+'.0d0/timeUnit    ',s)
 s=s.replace('ceiling(2000.0d0/dble(refineRatio))','ceiling(8.0d0/dble(refineRatio))')
 if steady:
  s=s.replace('!#define steadyFlow','#define steadyFlow').replace('#define unsteadyFlow','!#define unsteadyFlow')
  s=s.replace('ceiling(20000000.0d0/dble(refineRatio))','ceiling(38.0d0/dble(refineRatio))')
 if periodic:
  s=s.replace('#define VerticalWallsNoslip','!#define VerticalWallsNoslip').replace('!#define VerticalWallsPeriodicalU','#define VerticalWallsPeriodicalU')
  s=s.replace('#define VerticalWallsAdiabatic','!#define VerticalWallsAdiabatic').replace('!#define VerticalWallsPeriodicalT','#define VerticalWallsPeriodicalT')
 return s
results={}
for label,steady,periodic,ratio in [('wall',False,False,2),('periodic',False,True,2),('uniform_periodic',False,True,1),('steady',True,False,2)]:
 d=root/label;d.mkdir(exist_ok=True);(d/'solver.F90').write_text(small(src,steady,periodic,ratio),encoding='utf-8')
 c=subprocess.run([compiler,'-cpp','-ffree-line-length-none','-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow','solver.F90','-lgomp','-o','solver.exe'],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace")
 (d/'compile.log').write_text(c.stdout+c.stderr,encoding='utf-8');assert c.returncode==0,(label,c.stderr)
 r=subprocess.run([str(d/'solver.exe')],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=60)
 (d/'run.log').write_text(r.stdout+r.stderr,encoding='utf-8');results[label]={'code':r.returncode,'last':(r.stdout+r.stderr)[-1200:]}
assert all(v['code']==0 for v in results.values()),results
report={'smoke':results}

results={}
post=(repo/'多块网格/2DRBOpenaccMultiblockPostprocess.F90').read_text(encoding='utf-8')
def compile_run(d,s,stem='solver',expect=0):
 d.mkdir(exist_ok=True);(d/(stem+'.F90')).write_text(s,encoding='utf-8')
 c=subprocess.run([compiler,'-cpp','-ffree-line-length-none','-fcheck=all','-ffpe-trap=invalid,zero,overflow',stem+'.F90','-lgomp','-o',stem+'.exe'],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace")
 assert c.returncode==0,c.stderr
 r=subprocess.run([str(d/(stem+'.exe'))],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=60)
 (d/(stem+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8');assert r.returncode==expect,(str(d),r.stdout+r.stderr)
 return r
for label in ['wall','periodic','uniform_periodic','steady']:
 d=root/label;ps=post.replace('=500.0d0','=0.02d0').replace('=1000.0d0','=0.09d0')
 if label=='steady':ps=ps.replace('!#define steadyFlow','#define steadyFlow').replace('#define unsteadyFlow','!#define unsteadyFlow')
 compile_run(d,ps,'post')
 a=np.atleast_2d(np.loadtxt(d/'NuRe_2DOpenaccMultiblock.dat'));b=np.atleast_2d(np.loadtxt(d/'SnapshotNuRe_2DOpenaccMultiblock.dat'))
 assert a.shape[1]==15 and b.shape[1]==14
 assert np.all(a[:,5]<=a[:,4]+1e-12*np.maximum(1,np.abs(a[:,4])))
 match=b[np.isin(b[:,0],a[:,1])];diff=float(np.max(np.abs(a[:,2:]-match[:,1:])))
 np.testing.assert_allclose(a[:,2:],match[:,1:],rtol=1e-13,atol=1e-13)
 results[label+'_snapshot_max_diff']=diff
 if label=='steady':
  finaltext=(d/'FinalNuRe_2DOpenaccMultiblock.dat').read_text()
  assert '# itc_total = 38' in finaltext and '# error_check_itc = 32' in finaltext
  assert '# errorU errorT = ' in finaltext
  assert (d/'SteadySummary_2DOpenaccMultiblock.dat').read_text()==finaltext
 # Output index uses full filenames with the exact cumulative integer step.
 for row in (d/'OutputFiles_2DOpenaccMultiblock.dat').read_text().splitlines():
  step,name=row.split()
  if name.startswith('buoyancy'):
   assert name.endswith(f'tff-{step}itc{Path(name).suffix}')
   assert (d/name).is_file()
 # Snapshot 0002: full arrays/original offsets; 11 ints, 8 reals, box(4), block geometry, 4 fields.
 meta=(d/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().splitlines()
 settings=list(map(int,meta[3].split()))
 meta_step=int(meta[meta.index('# itc_total NuReSampleCount t_ff errorU errorT')+1].split()[0])
 expected_geometry=[]
 for block in (['coarse','left','right','bottom','top'] if settings[-1]==5 else ['coarse']):
  pos=next(i for i,line in enumerate(meta) if line.startswith('# '+block+':'))
  dims=list(map(int,meta[pos+1].split()))
  geometry=list(map(float,meta[pos+2].split()))
  expected_geometry.append((dims[:2],geometry[:3]))
 for snap in d.glob('buoyancy*Snapshot-*.bin'):
  with snap.open('rb') as f:
   assert f.read(16)==b'MB2DSNAPSHOT0002'
   nb,nx,ny,total_step,mode,ratio,skin,left,right,bottom,top=struct.unpack('11i',f.read(44))
   assert [nx,ny,ratio,left,right,bottom,top,skin,nb]==settings
   physics=struct.unpack('8d',f.read(64));box=struct.unpack('4d',f.read(32));cells=0;block_cells=[]
   assert abs(physics[0]-total_step/physics[7])<1e-12*max(1,abs(physics[0]))
   for i in range(nb):
    ni,nj=struct.unpack('2i',f.read(8));offsets=struct.unpack('3d',f.read(24))
    assert [ni,nj]==expected_geometry[i][0]
    np.testing.assert_array_equal(offsets,expected_geometry[i][1])
    block_cells.append(ni*nj);cells+=ni*nj
   fields=np.frombuffer(f.read(),dtype=np.float64).reshape(4,cells)
  assert snap.stat().st_size==16+44+96+32*nb+32*cells
  if total_step==meta_step:
   # The full coarse u/v/T/rho, including receivers, must equal the checkpoint arrays exactly.
   checkpoint_fields=np.frombuffer((d/meta[1]).read_bytes(),dtype=np.float64,
       count=4*block_cells[0],offset=14*block_cells[0]*8).reshape(4,block_cells[0])
   np.testing.assert_array_equal(fields[:,:block_cells[0]],checkpoint_fields)
 results[label+'_snapshot_layout']='full dimensions and original offsets match meta; coarse fields match checkpoint exactly'
for label in ['wall','periodic','steady']:
 base=root/label;d=root/(label+'_restart');shutil.copytree(base,d,dirs_exist_ok=True)
 latest='reloadFile2DOpenaccMultiblock-latest.meta'
 assert [p.name for p in base.glob('reload*.meta')]==[latest]
 # Stop a separate run after its first checkpoint; only latest.meta is produced.
 checkpoint=root/(label+'_checkpoint22')
 first=(base/'solver.F90').read_text(encoding='utf-8').replace('  end subroutine output_ReloadFile',
     "    if(reloadItcOffset+itc==22) stop\n  end subroutine output_ReloadFile")
 compile_run(checkpoint,first)
 meta22=(checkpoint/latest).read_bytes()
 chosen=meta22.decode().splitlines()[1]
 assert chosen.endswith('tff-22itc.bin')
 (d/latest).write_bytes(meta22)
 # Existing latest must override even an invalid manual filename/step.
 s=(base/'solver.F90').read_text(encoding='utf-8').replace('loadInitField=0','loadInitField=1').replace('reloadFileName=""','reloadFileName="ignored.bin"')
 s=s.replace(':: reloadItcOffset=0',':: reloadItcOffset=-1',1)
 if label=='steady':s=s.replace('ceiling(38.0d0/dble(refineRatio))','ceiling(16.0d0/dble(refineRatio))')
 if label=='wall':
  initOnly=s.replace('    call initial()','    call initial()\n    stop',1)
  compile_run(d,initOnly)
  assert (d/latest).read_bytes()==meta22
  results['latest_priority']='latest metadata overrides manual filename/step; initialization preserves it'
 # A damaged discarded tail must not block a valid checkpoint restart.
 with (d/'NuRe_2DOpenaccMultiblock.dat').open('a') as f:f.write('incomplete discarded tail\n')
 compile_run(d,s)
 final=(d/latest).read_text().splitlines()[1]
 expected_step=38 if label=='steady' else 60
 assert final.endswith(f'tff-{expected_step}itc.bin')
 ref=(base/latest).read_text().splitlines()[1]
 assert (d/final).read_bytes()==(base/ref).read_bytes(),label+' restart arrays differ'
 assert (d/'NuRe_2DOpenaccMultiblock.dat').read_bytes()==(base/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
 if label=='steady':np.testing.assert_array_equal(np.loadtxt(d/'Convergence_2DOpenaccMultiblock.dat'),np.loadtxt(base/'Convergence_2DOpenaccMultiblock.dat'))
 results[label+'_restart']='bitwise identical arrays and history'
 # Missing latest: the same binary plus manually supplied integer step restores the state.
 manual=root/(label+'_manual');shutil.copytree(base,manual);(manual/latest).unlink()
 ms=s.replace('reloadFileName="ignored.bin"',f'reloadFileName="{chosen}"').replace(':: reloadItcOffset=-1',':: reloadItcOffset=22',1)
 compile_run(manual,ms)
 assert (manual/ref).read_bytes()==(base/ref).read_bytes()
 assert (manual/'NuRe_2DOpenaccMultiblock.dat').read_bytes()==(base/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
 assert [p.name for p in manual.glob('reload*.meta')]==[latest]
 results[label+'_manual_restart']='missing latest: bitwise identical arrays/history; only latest.meta saved'
 if label=='wall':
  # Existing invalid metadata must never fall back to otherwise valid manual settings.
  for kind in ['corrupt','truncated','physics_mismatch','negative_manual_step','unsynchronized_manual_step','missing_manual_file']:
   bad=root/kind;shutil.copytree(base,bad);bs=ms
   if kind=='corrupt':(bad/latest).write_text('broken metadata\n')
   elif kind=='truncated':(bad/latest).write_text('MB2D_RESTART 1\n')
   elif kind=='physics_mismatch':
    rows=meta22.decode().splitlines();pi=next(i for i,r in enumerate(rows) if r.startswith('# Rayleigh'))+1
    vals=rows[pi].split();vals[0]='1.0';rows[pi]=' '.join(vals);(bad/latest).write_text('\n'.join(rows)+'\n')
   else:
    (bad/latest).unlink()
    if kind=='negative_manual_step':bs=ms.replace(':: reloadItcOffset=22',':: reloadItcOffset=-2',1)
    elif kind=='unsynchronized_manual_step':bs=ms.replace(':: reloadItcOffset=22',':: reloadItcOffset=21',1)
    else:bs=ms.replace(f'reloadFileName="{chosen}"','reloadFileName="missing.bin"')
   original_history=(bad/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
   compile_run(bad,bs,expect=1)
   assert (bad/'NuRe_2DOpenaccMultiblock.dat').read_bytes()==original_history
   results[kind]='rejected before history changes'
  # No new steps: existing snapshot/plt names must be found, without duplicate index rows.
  zero=root/'zero_advance';shutil.copytree(base,zero)
  original_index=(zero/'OutputFiles_2DOpenaccMultiblock.dat').read_bytes()
  compile_run(zero,s)
  assert (zero/'OutputFiles_2DOpenaccMultiblock.dat').read_bytes()==original_index
  results['zero_advance_restart']='existing snapshot/plt found; no duplicate output rows'
  # Restore earlier latest, end at 36, then extend to 84 using the newly saved latest.
  (d/latest).write_bytes(meta22)
  s36=s.replace('unsteadyRunDuration=60.0d0/timeUnit','unsteadyRunDuration=36.0d0/timeUnit')
  stale=[p for p in d.glob('buoyancy*') if p.is_file() and p.name.endswith('-60itc'+p.suffix)]
  compile_run(d,s36)
  index=(d/'OutputFiles_2DOpenaccMultiblock.dat').read_text().splitlines()
  assert all(int(row.split()[0])<=36 for row in index if row)
  assert all(not p.exists() for p in stale),stale
  s84=s.replace('unsteadyRunDuration=60.0d0/timeUnit','unsteadyRunDuration=84.0d0/timeUnit')
  compile_run(d,s84)
  fresh=root/'wall84';freshsrc=(base/'solver.F90').read_text(encoding='utf-8').replace('unsteadyRunDuration=60.0d0/timeUnit','unsteadyRunDuration=84.0d0/timeUnit')
  compile_run(fresh,freshsrc)
  final=(d/latest).read_text().splitlines()[1]
  assert (d/final).read_bytes()==(fresh/final).read_bytes()
  assert (d/'NuRe_2DOpenaccMultiblock.dat').read_bytes()==(fresh/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
  results['extend_target_and_delete_future']='passed'
report['restart_and_snapshot']=results

if args.before:
    # The added diagnostic must leave fields and every pre-existing history column unchanged.
    for label,steady,periodic,ratio in [('wall',False,False,2),('periodic',False,True,2),('uniform_periodic',False,True,1),('steady',True,False,2)]:
        d=root/('original_'+label)
        compile_run(d,small(args.before.read_text(encoding='utf-8'),steady,periodic,ratio))
        fresh=root/label
        fn=(d/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().splitlines()[1]
        fn2=(fresh/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().splitlines()[1]
        assert (d/fn).read_bytes()==(fresh/fn2).read_bytes(),'Pre-edit numerical state changed'
        old=np.atleast_2d(np.loadtxt(d/'NuRe_2DOpenaccMultiblock.dat'))
        new=np.atleast_2d(np.loadtxt(fresh/'NuRe_2DOpenaccMultiblock.dat'))
        np.testing.assert_array_equal(old,np.delete(new,5,axis=1))
    report['before_edit_comparison']='four cases: bitwise identical checkpoint fields and all original Nu/Re/history columns'

# Independent area-weighted diagnostic: coarse speed=5, fine speed=1.
# The physical center area gives analytic mean and RMS; uniform mode gives speed=5.
for ratio in [1,2]:
    s=small(src,steady=True,ratio=ratio)
    s=re.sub(r'    program main.*?    end program main','',s,flags=re.S)
    probe='''
program mean_re_probe
 use commondata
 implicit none
 real(kind=8) :: nu,re,reMean,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax,fraction,expected
 call initial()
 u_coarse=3d0
 v_coarse=4d0
 fraction=1d0
 if(refineRatio>1) then
   u_left=0.6d0; v_left=0.8d0
   u_right=0.6d0; v_right=0.8d0
   u_bottom=0.6d0; v_bottom=0.8d0
   u_top=0.6d0; v_top=0.8d0
   fraction=(centerIntegrationBox(2)-centerIntegrationBox(1))* &
       (centerIntegrationBox(4)-centerIntegrationBox(3))/(dble(nx)*dble(ny))
 endif
 itc=outputSnapshotIntervalItc
 call calculate_NuRe(nu,re,reMean,hot,cold,middle,mass,meanT,tmin,tmax,rmin,rmax)
 expected=(1d0+4d0*fraction)*lengthUnit/viscosity
 if(abs(reMean-expected)>1d-12*max(1d0,abs(expected))) error stop 'Mean Re area weighting'
 expected=sqrt(1d0+24d0*fraction)*lengthUnit/viscosity
 if(abs(re-expected)>1d-12*max(1d0,abs(expected))) error stop 'RMS Re area weighting'
 call calNuRe()
 call output_SnapshotFile()
 call output_final_NuRe()
end program mean_re_probe
'''
    d=root/f'mean_re_ratio{ratio}'
    compile_run(d,s+probe)
    compile_run(d,post.replace('!#define steadyFlow','#define steadyFlow').replace('#define unsteadyFlow','!#define unsteadyFlow'),'post')
    history=np.atleast_2d(np.loadtxt(d/'NuRe_2DOpenaccMultiblock.dat'))
    snapshot=np.atleast_2d(np.loadtxt(d/'SnapshotNuRe_2DOpenaccMultiblock.dat'))
    final=np.atleast_2d(np.loadtxt(d/'FinalNuRe_2DOpenaccMultiblock.dat'))
    np.testing.assert_allclose(history[:,2:],snapshot[:,1:],rtol=1e-13,atol=1e-13)
    np.testing.assert_array_equal(history[:,2:],final)
report['mean_re_analytic']='coarse/fine physical-area weighting and uniform constant-speed limit; history, snapshot and final results agree'

# Standalone halo probe with unique direction/row/column markers.
s=(root/'periodic'/'solver.F90').read_text(encoding='utf-8');s=re.sub(r'    program main.*?    end program main', '',s,flags=re.S)
probe='''
program periodic_probe
 use commondata
 implicit none
 integer :: i,j,a,srcj,nr
 real(kind=8) :: expect
 call initial()
 nr=size(f_right,1)
 do a=0,8
   do j=1,ny
     do i=1,size(f_left,1)
       f_post_left(i,j,a)=10000d0+100d0*a+2d0*j+0.01d0*i
     enddo
     do i=1,nr
       f_post_right(i,j,a)=20000d0+100d0*a+2d0*j+0.01d0*i
     enddo
   enddo
 enddo
 call exchange_fine_flow_halos()
 call streaming(size(f_left,1),ny,f_left,f_post_left)
 call streaming(nr,ny,f_right,f_post_right)
 call bounceback(size(f_left,1),ny,f_left,f_post_left,.true.,.false.,.true.,.true.)
 call bounceback(nr,ny,f_right,f_post_right,.false.,.true.,.true.,.true.)
 do j=1,ny
   do a=1,8
     srcj=j-ey(a)
     if(srcj<1 .OR. srcj>ny) cycle ! Those links hit the horizontal wall.
     if(ex(a)==1) then
       expect=20000d0+100d0*a+2d0*srcj+0.01d0*nr
       if(f_left(1,j,a)/=expect) error stop 'Wrong left periodic link'
     else if(ex(a)==-1) then
       expect=10000d0+100d0*a+2d0*srcj+0.01d0
       if(f_right(nr,j,a)/=expect) error stop 'Wrong right periodic link'
     endif
   enddo
 enddo
 if(any(f_left(:,1,5)/=f_post_left(1:size(f_left,1),1,7))) error stop 'Bottom corner reflection'
 if(any(f_right(:,ny,7)/=f_post_right(1:nr,ny,5))) error stop 'Top corner reflection'
 do a=0,4
   do j=1,ny
     g_post_left(1,j,a)=100d0*a+j
     g_post_right(nr,j,a)=200d0*a+j
   enddo
 enddo
 call exchange_fine_thermal_halos()
 call streamingT(size(g_left,1),ny,g_left,g_post_left)
 call streamingT(nr,ny,g_right,g_post_right)
 call bouncebackT(size(g_left,1),ny,g_left,g_post_left,.true.,.false.,.true.,.true.)
 call bouncebackT(nr,ny,g_right,g_post_right,.false.,.true.,.true.,.true.)
 do j=1,ny
   if(g_left(1,j,1)/=200d0+j) error stop 'Wrong left thermal link'
   if(g_right(nr,j,3)/=300d0+j) error stop 'Wrong right thermal link'
 enddo
 print *, 'Periodic directions, vertical shifts and wall corners passed.'
end program periodic_probe
'''
d=root/'haloprobe';d.mkdir(exist_ok=True);(d/'probe.F90').write_text(s+probe,encoding='utf-8')
c=subprocess.run([fc,'-cpp','-ffree-line-length-none','-fcheck=all','probe.F90','-lgomp','-o','probe.exe'],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace");assert c.returncode==0,c.stderr
r=subprocess.run([str(d/'probe.exe')],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=60);assert r.returncode==0,r.stderr
report['periodic_links']=r.stdout.strip()

base=(root/'wall'/'solver.F90').read_text(encoding='utf-8')
post=Path('多块网格/2DRBOpenaccMultiblockPostprocess.F90').read_text(encoding='utf-8').replace('=500.0d0','=0.02d0').replace('=1000.0d0','=0.09d0')
results={}
variants={
 'velocity_periodic':base.replace('#define VerticalWallsNoslip','!#define VerticalWallsNoslip').replace('!#define VerticalWallsPeriodicalU','#define VerticalWallsPeriodicalU'),
 'thermal_periodic':base.replace('#define VerticalWallsAdiabatic','!#define VerticalWallsAdiabatic').replace('!#define VerticalWallsPeriodicalT','#define VerticalWallsPeriodicalT'),
 'ratio4':base.replace('nx=32, ny=32','nx=64, ny=64').replace('refineRatio=2','refineRatio=4').replace('fineLayerCellsLeft=8','fineLayerCellsLeft=16').replace('fineLayerCellsRight=9','fineLayerCellsRight=17').replace('fineLayerCellsBottom=8','fineLayerCellsBottom=16').replace('fineLayerCellsTop=9','fineLayerCellsTop=17'),
 'sideheated':base.replace('#define RayleighBenardCell','!#define RayleighBenardCell').replace('!#define SideHeatedCell','#define SideHeatedCell').replace('#define HorizontalWallsConstT','!#define HorizontalWallsConstT').replace('!#define HorizontalWallsAdiabatic','#define HorizontalWallsAdiabatic').replace('#define VerticalWallsAdiabatic','!#define VerticalWallsAdiabatic').replace('!#define VerticalWallsConstT','#define VerticalWallsConstT')}
for label,s in variants.items():
 d=root/label;d.mkdir(exist_ok=True);(d/'solver.F90').write_text(s,encoding='utf-8');(d/'post.F90').write_text(post,encoding='utf-8')
 for stem in ['solver','post']:
  c=subprocess.run([fc,'-cpp','-ffree-line-length-none','-fcheck=all','-ffpe-trap=invalid,zero,overflow',stem+'.F90','-lgomp','-o',stem+'.exe'],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace");assert c.returncode==0,c.stderr
  r=subprocess.run([str(d/(stem+'.exe'))],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=60);assert r.returncode==0,r.stderr
 import numpy as np
 a=np.loadtxt(d/'NuRe_2DOpenaccMultiblock.dat');b=np.loadtxt(d/'SnapshotNuRe_2DOpenaccMultiblock.dat');b=b[np.isin(b[:,0],a[:,1])]
 np.testing.assert_allclose(a[:,2:],b[:,1:],rtol=1e-12,atol=1e-12)
 results[label]='run and snapshot diagnostics passed'
# Compile all actual OpenACC directives, separately from the CPU probes.
d=root/'accsyntax';d.mkdir(exist_ok=True);(d/'solver.F90').write_text(Path('多块网格/2DRBOpenaccMultiblock.F90').read_text(encoding='utf-8'),encoding='utf-8')
c=subprocess.run([fc,'-cpp','-ffree-line-length-none','-fopenacc','-fsyntax-only','solver.F90'],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace")
(d/'compile.log').write_text(c.stdout+c.stderr,encoding='utf-8');results['gfortran_openacc_syntax']={'exit':c.returncode,'log':c.stderr[:1400]}
assert c.returncode==0,c.stderr
report['variants']=results

# Independent analytic statistics check (linear Nu and linear Re^2).
d=root/'analytic';d.mkdir()
ps=(repo/'多块网格/2DRBOpenaccMultiblockPostprocess.F90').read_text(encoding='utf-8')
ps=ps.replace('processSnapshots=.true.','processSnapshots=.false.').replace('=500.0d0','=0.5d0').replace('=1000.0d0','=2.5d0')
rows=[]
for i in range(1,6):
 t=0.5*i
 # Deliberately distinct ReVolMean: window statistics must skip this column.
 rows.append(f'{i} {10*i} {t} {2*t} {(1+t)**0.5} {0.1*t} 7 8 9 1 0 -0.5 0.5 1 1')
(d/'NuRe_2DOpenaccMultiblock.dat').write_text('\n'.join(rows)+'\n')
compile_run(d,ps,'post')
line=next(x for x in (d/'NuReStatistics_2DOpenaccMultiblock.dat').read_text().splitlines() if x.startswith('whole:'))
np.testing.assert_allclose(np.fromstring(line[6:],sep=' '),[3,2.5**0.5,7,8,9],rtol=1e-14,atol=1e-14)
(d/'NuRe_2DOpenaccMultiblock.dat').write_text('\n'.join(rows[:2]+rows[3:])+'\n')
compile_run(d,ps,'post',expect=1)
report['analytic_window']='linear Nu and RMS Re verified; ReVolMean excluded from five-column summary; missing sample rejected'

# A total target between samples rounds upward; a positive sub-coarse interval is accepted.
d=root/'aligned_target'
s=(root/'wall'/'solver.F90').read_text(encoding='utf-8').replace('unsteadyRunDuration=60.0d0/timeUnit','unsteadyRunDuration=61.0d0/timeUnit')
compile_run(d,s)
meta=(d/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().splitlines()
assert int(meta[meta.index('# itc_total NuReSampleCount t_ff errorU errorT')+1].split()[0])==66
smallInterval=s.replace('outputSnapshotInterval=6.0d0/timeUnit','outputSnapshotInterval=0.1d0/timeUnit').replace('unsteadyRunDuration=61.0d0/timeUnit','unsteadyRunDuration=3.0d0/timeUnit')
d=root/'subcoarse_interval';compile_run(d,smallInterval)
a=np.atleast_2d(np.loadtxt(d/'NuRe_2DOpenaccMultiblock.dat'));np.testing.assert_array_equal(a[:,1],[2,4])
report['alignment']='61 requested steps -> 66; positive 0.1 step interval -> 2'

# Rounded display times may coincide; the integer step still makes filenames unique.
s=(root/'wall'/'solver.F90').read_text(encoding='utf-8').replace('timeUnit=sqrt(lengthUnit/gBeta)','timeUnit=1.0d14')
for key,old in [('outputSnapshotInterval',6),('outputPltFileInterval',12),('reloadFileInterval',20)]:
 s=s.replace(f'{key}={old}.0d0/timeUnit',f'{key}=2.0d0/timeUnit')
s=s.replace('unsteadyRunDuration=60.0d0/timeUnit','unsteadyRunDuration=6.0d0/timeUnit')
d=root/'rounded_display_time';compile_run(d,s)
for pattern in ['buoyancy*Snapshot-*.bin','buoyancy*Tecplot-*.dat','reload*-*itc.bin']:
 files=list(d.glob(pattern));assert len(files)>=2
 assert all('-0tff-' in p.name for p in files)
 assert len({p.name for p in files})==len(files)
report['rounded_display_time']='identical displayed tff retains distinct snapshots, plt and checkpoints by step'


# Invalid checkpoint state must be rejected before the history is truncated.
d=root/'invalid_restart';shutil.copytree(root/'wall',d)
metaPath=d/'reloadFile2DOpenaccMultiblock-latest.meta'
lines=metaPath.read_text().splitlines();i=lines.index('# itc_total NuReSampleCount t_ff errorU errorT')+1
parts=lines[i].split();parts[0]='59';lines[i]=' '.join(parts);metaPath.write_text('\n'.join(lines)+'\n')
history=(d/'NuRe_2DOpenaccMultiblock.dat').read_bytes()
compile_run(d,(root/'wall'/'solver.F90').read_text(encoding='utf-8').replace('loadInitField=0','loadInitField=1'),expect=1)
assert (d/'NuRe_2DOpenaccMultiblock.dat').read_bytes()==history
report['invalid_checkpoint']='unsynchronized step rejected without rewriting history'

# Actual directives executed on libgomp host; not a NVIDIA GPU benchmark.
d=root/'acc_host';d.mkdir()
(d/'solver.F90').write_text((root/'periodic'/'solver.F90').read_text(encoding='utf-8'),encoding='utf-8')
c=subprocess.run([fc,'-cpp','-ffree-line-length-none','-fopenacc','-fcheck=all','solver.F90','-o','solver.exe'],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace")
assert c.returncode==0,c.stderr
r=subprocess.run([str(d/'solver.exe')],cwd=d,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=60)
assert r.returncode==0,r.stderr
ref=root/'periodic';fname=(d/'reloadFile2DOpenaccMultiblock-latest.meta').read_text().splitlines()[1]
assert (d/fname).read_bytes()==(ref/fname).read_bytes()
report['openacc_host']='directives executed; state identical to CPU reference'

(root/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
