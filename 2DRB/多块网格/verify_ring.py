"""Connected-ring topology, transport, thermal and restart regression on OpenACC host.
Temporary source overrides only; this script never edits the solver.
"""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import verify_multiblock as v

GEOMETRY = r"""program main
 use commondata
 use openacc
 implicit none
 integer :: b,i,j,l,ii,jj
 real(8) :: area,a
 call acc_init(acc_device_host)
 call initial()
 if(nBlocks/=2 .or. nLinks/=2) error stop 'Expected only two cross-level links'
 do l=1,nLinks
  if(linkReceiver(l)==linkDonor(l)) error stop 'Self link'
  if(blockH(linkReceiver(l))==blockH(linkDonor(l))) error stop 'Same-level exchange remains'
 enddo
 area=0d0
 do b=1,nBlocks
  do j=1,blockNj(b)
   do i=1,blockNi(b)
    a=owned_cell_area(b,i,j)
    if(a<0d0) error stop 'Negative ownership'
    if(b==2 .and. a>0d0) then
     if(.not.fine_active(i,j)) error stop 'Inactive owned node'
    endif
    area=area+a
   enddo
  enddo
 enddo
 if(abs(area-dble(nx)*ny)>1d-10) error stop 'Area mismatch'
 ! Direct migration across the former upper-left fine seam (no exchange call).
 call select_block(2)
 f_post=0d0;g_post=0d0
 f(nx/2,ny/2,0)=0.314159d0;g(nx/2,ny/2,0)=0.271828d0
 ii=4;jj=ny-fineLayerCellsTop+1
 f_post(ii,jj,2)=0.123456789d0
 g_post(ii,jj,2)=0.987654321d0
 call enter_data_2d_openacc()
 ! enter_data initializes post arrays on device; now explicitly load the marker fields.
 !$acc update device(f_post,g_post)
 call streaming(nx,ny,f,f_post)
 call streamingT(nx,ny,g,g_post)
 call update_host_all(.true.)
 call select_block(2)
 if(f(ii,jj+1,2)/=0.123456789d0) error stop 'Flow seam migration failed'
 if(g(ii,jj+1,2)/=0.987654321d0) error stop 'Thermal seam migration failed'
 if(f(nx/2,ny/2,0)/=0.314159d0 .or. g(nx/2,ny/2,0)/=0.271828d0) error stop 'Inactive hole advanced'
 call exit_data_2d_openacc()
 write(*,*) 'RING_GEOMETRY_AND_DIRECT_STREAMING_OK',area
end program main"""


def main():
 s=v.SOURCE.read_text(encoding='utf-8-sig')
 report={'source_sha256':hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
         'device':'gfortran OpenACC host','build_directory':str(v.BUILD),'checks':[]}
 for name,kw in [('rb',{}),('legacy',{'legacy':True}),('side',{'side':True}),('steady',{'steady':True}),('ha',{'side':True,'ha':True})]:
  v.compile_source('syntax_'+name,v.variant(s,**kw),syntax=True)
 report['checks'].append({'five_macro_syntax_checks':True})
 for ratio in [2,4,8]:
  folder,exe=v.compile_source('geometry_'+str(ratio),v.variant(s,ratio=ratio),GEOMETRY,n=192,ny=160,walls=(32,33,32,33),ra=10000)
  output=v.run([exe],folder)
  report['checks'].append({'ratio':ratio,'geometry_and_direct_seam_streaming':output.strip()})
  print('Geometry and direct streaming',ratio,flush=True)
 for name,kw in [('rb',{}),('side',{'side':True}),('legacy',{'side':True,'legacy':True})]:
  d=v.CONDUCTION.replace('                err=max(err','                if(owned_cell_area(b,i,j)<=0d0) cycle\n                err=max(err')
  folder,exe=v.compile_source('conduction_'+name,v.variant(s,**kw),d,ra=10000)
  output=v.run([exe],folder,timeout=150)
  values=[float(x) for x in next(l for l in output.splitlines() if 'CONDUCTION_ERROR' in l).split()[1:]]
  hist=np.loadtxt(folder/'NuRe_2DOpenaccMultiblock.dat')
  assert np.max(abs(hist[[1,3,4,5]]-1))<1e-10
  report['checks'].append({'conduction':name,'coarse_steps':300,'max_T_and_rho_error':values,'Nu':hist[[1,3,4,5]].tolist()})
  print('Conduction',name,values,flush=True)
 for name,kw in [('rb',{}),('legacy',{'side':True,'legacy':True})]:
  parent=v.PARENT.read_text(encoding='utf-8-sig')
  a,ae=v.compile_source('uniform_parent_'+name,v.variant(parent,**kw),v.uniform_driver(True),n=32)
  b,be=v.compile_source('uniform_ring_'+name,v.variant(s,ratio=1,**kw),v.uniform_driver(),n=32)
  v.run([ae],a);v.run([be],b)
  diff=np.max(abs(np.fromfile(a/'state.bin',dtype='<f8')-np.fromfile(b/'state.bin',dtype='<f8')))
  assert diff<3e-13
  report['checks'].append({'uniform_regression':name,'max_state_error':float(diff)})
 for name,kw in [('rb',{}),('steady',{'steady':True})]:
  full,fe=v.compile_source('restart_full_'+name,v.variant(s,**kw),v.RESTART_DRIVER)
  split,se=v.compile_source('restart_split_'+name,v.variant(s,**kw),v.RESTART_DRIVER)
  _,re=v.compile_source('restart_resume_'+name,v.variant(s,restart=True,**kw),v.RESTART_DRIVER)
  v.run([fe,'40'],full);v.run([se,'20'],split);v.run([re,'40'],split)
  a=np.fromfile(full/'allstate.bin',dtype='<f8');b=np.fromfile(split/'allstate.bin',dtype='<f8')
  assert np.all(np.isfinite(a)) and np.array_equal(a,b)
  report['checks'].append({'exact_restart':name,'fine_steps':40,'split':20})
  print('Exact restart',name,flush=True)
 report['status']='passed'
 (v.HERE/'ring_verification_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print('All ring checks passed',flush=True)

if __name__=='__main__': main()

