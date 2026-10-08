"""Verify one reusable transfer buffer and separate per-region time history.

Usage: python verify_shared_buffer.py PATH_TO_PRE_CHANGE_SOURCE
All builds and numerical output stay in a temporary directory (OpenACC host).
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import verify_multiblock as v
import verify_simplified as s
import verify_direct_exchange as d
import verify_packed_buffer as b

def driver(shared, steps):
    body=b.driver(steps)
    if shared:
        body=body.replace('integer :: step,k','integer :: step,k,region,iStart,iEnd,jStart,jEnd\n    real(8) :: xOffset,yOffset')
        body=body.replace('    write(78) coarseBuffer',"""    if(refineRatio==1) then
        write(78) coarseHistory
    else
        do k=0,2
            do region=1,4
                call coarse_receiver_layout(region,iStart,iEnd,jStart,jEnd,xOffset,yOffset)
                write(78) coarseHistory(1:20*(iEnd-iStart+1)*(jEnd-jStart+1),k,region)
            enddo
        enddo
    endif""")
    return body

def reuse_driver():
    body="""program main
    use commondata
    use openacc
    implicit none
    integer :: p,k,region,i,j,iStart,iEnd,jStart,jEnd,n
    real(8) :: xOffset,yOffset,expected
    call acc_init(acc_device_host)
    call initial()
    call enter_data_2d_openacc()
    !$acc wait(1)
    do region=1,4
        do k=0,2
            do p=1,size(coarseHistory,1)
                coarseHistory(p,k,region)=0.01d0*region+0.001d0*k+0.0000001d0*p
            enddo
        enddo
    enddo
    coarseBuffer=-999.0d0
    open(unit=78,file='history_before.bin',form='unformatted',access='stream',status='replace')
    write(78) coarseHistory
    close(78)
    !$acc update device(coarseHistory,coarseBuffer) async(1)
"""
    for name in d.REGIONS[1:]:
        body+=f"""    rho_{name}=-123.0d0
    T_{name}=-456.0d0
    !$acc update device(rho_{name},T_{name}) async(1)
"""
    for region,name in enumerate(d.REGIONS[1:],1):
        body+=f"""    call interpolate_coarse_to_fine({region},size(f_{name},1),size(f_{name},2), &
        f_{name},g_{name},rho_{name},u_{name},v_{name},T_{name},Fx_{name},Fy_{name},[0.25d0,0.5d0,0.25d0])
    call coarse_receiver_layout({region},iStart,iEnd,jStart,jEnd,xOffset,yOffset)
    n=20*(iEnd-iStart+1)*(jEnd-jStart+1)
    !$acc update self(coarseBuffer) async(1)
    !$acc wait(1)
    do p=1,size(coarseBuffer)
        expected=0.0d0
        if(p<=n) expected=0.01d0*{region}+0.001d0+0.0000001d0*p
        if(abs(coarseBuffer(p)-expected)>1.0d-14) error stop 'Shared buffer content or cleared tail differs'
    enddo
"""
    body+="""    call update_host_all(.true.)
    open(unit=78,file='history_after.bin',form='unformatted',access='stream',status='replace')
    write(78) coarseHistory
    close(78)
"""
    for region,name in enumerate(d.REGIONS[1:],1):
        body+=f"""    call coarse_receiver_layout({region},iStart,iEnd,jStart,jEnd,xOffset,yOffset)
    do j=1,size(f_{name},2)
        do i=1,size(f_{name},1)
            if(i>=iStart .AND. i<=iEnd .AND. j>=jStart .AND. j<=jEnd) then
                p=20*((j-jStart)*(iEnd-iStart+1)+i-iStart)
                expected=0.01d0*{region}+0.001d0+0.0000001d0*(p+1)
                if(abs(rho_{name}(i,j)-expected)>1.0d-14) error stop 'Wrong region or node received'
                expected=0.01d0*{region}+0.001d0+0.0000001d0*(p+4)
                if(abs(T_{name}(i,j)-expected)>1.0d-14) error stop 'Wrong thermal component received'
            else
                if(rho_{name}(i,j)/=-123.0d0 .OR. T_{name}(i,j)/=-456.0d0) &
                    error stop 'Non-receiver changed'
            endif
        enddo
    enddo
"""
    body+="""    call exit_data_2d_openacc()
    call deallocate_grid_arrays()
    if(allocated(coarseHistory) .OR. allocated(coarseBuffer)) error stop 'Exchange storage not released'
end program main"""
    return body

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline',type=Path)
    args=parser.parse_args()
    old=args.baseline.read_text(encoding='utf-8-sig')
    new=v.SOURCE.read_text(encoding='utf-8-sig')
    assert 'coarseBuffer(:,:)' in old and 'coarseHistory(:,:,:)' in new
    report=dict(source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest(),
                baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
                baseline_path=str(args.baseline.resolve()),device='gfortran OpenACC host',
                build_directory=str(v.BUILD),checks=[])
    cases=[
        ('uniform_rb',1,False,False,False,(16,17,16,17)),
        ('r2_rb',2,False,False,False,(16,17,16,17)),
        ('r2_side_steady',2,True,True,False,(16,17,16,17)),
        ('r4_side',4,True,False,False,(16,17,16,17)),
        ('r8_rb_steady',8,False,True,False,(32,33,32,33)),
        ('r2_side_ha',2,True,True,True,(16,17,16,17)),
        ('r2_side_seam',2,True,False,False,(44,13,12,13)),
    ]
    try:
        for label,ratio,side,steady,ha,walls in cases:
            dims=dict(n=128,ny=112,ra=10000) if ratio==8 else dict(n=96,ny=80,ra=10000)
            dims['walls']=walls
            settings=dict(ratio=ratio,side=side,steady=steady,ha=ha)
            folders=[]
            for tag,source,shared in [('before',old,False),('after',new,True)]:
                folder,exe=v.compile_source(label+'_'+tag,s.short_settings(source,**settings),driver(shared,40),**dims)
                v.run([exe],folder);folders.append(folder)
            before,after=folders
            s.compare_values(d.canonical(before),d.canonical(after),label+' fields',exact=True)
            s.compare_values(b.packets(before),b.packets(after),label+' history',exact=True)
            diagnostics=['NuRe_2DOpenaccMultiblock.dat']
            if steady: diagnostics+=['Convergence_2DOpenaccMultiblock.dat']
            for name in diagnostics:
                assert (before/name).read_bytes()==(after/name).read_bytes(),(label,name)
            s.compare_snapshots(next(before.glob('*Snapshot-*.bin')),next(after.glob('*Snapshot-*.bin')))
            assert d.c.checkpoint(before).read_bytes()==d.c.checkpoint(after).read_bytes(),label+' checkpoint bytes'
            restored,exe=v.compile_source(label+'_restore',s.short_settings(new,restart=True,**settings),driver(True,0),**dims)
            d.copy_checkpoint(before,restored);v.run([exe],restored)
            s.compare_values(d.canonical(before),d.canonical(restored),label+' restore',exact=True)
            s.compare_values(b.packets(before),b.packets(restored),label+' restore history',exact=True)
            half=(20//ratio)*ratio
            for tag,first,first_shared,second,second_shared in [
                ('old_to_new',old,False,new,True),
                ('new_to_old',new,True,old,False),
                ('new_to_new',new,True,new,True)]:
                split,exe=v.compile_source(label+'_'+tag+'_half',s.short_settings(first,**settings),driver(first_shared,half),**dims)
                v.run([exe],split)
                _,exe=v.compile_source(label+'_'+tag+'_resume',s.short_settings(second,restart=True,**settings),driver(second_shared,40-half),**dims)
                v.run([exe],split)
                s.compare_values(d.canonical(after),d.canonical(split),label+' '+tag,exact=True)
                s.compare_values(b.packets(after),b.packets(split),label+' '+tag+' history',exact=True)
                assert d.c.checkpoint(after).read_bytes()==d.c.checkpoint(split).read_bytes()
            entry=dict(case=label,fields_exact=True,all_history_components_exact=True,diagnostics_exact=True,
                       checkpoint_bytes_exact=True,old_checkpoint_restore_exact=True,
                       old_to_new_resume_exact=True,new_to_old_resume_exact=True,new_to_new_resume_exact=True)
            if ratio>1:
                poly,exe=v.compile_source(label+'_polynomial',s.short_settings(new,**settings),
                                         b.polynomial_driver().replace('save_coarse_buffer','save_coarse_history'),**dims)
                v.run([exe],poly)
                count,error=np.loadtxt(poly/'polynomial.dat')
                entry.update(polynomial_receivers=int(count),polynomial_max_abs_error=float(error))
                reuse,exe=v.compile_source(label+'_reuse',s.short_settings(new,**settings),reuse_driver(),**dims)
                v.run([exe],reuse)
                assert (reuse/'history_before.bin').read_bytes()==(reuse/'history_after.bin').read_bytes()
                entry.update(reuse_each_region_from_index_one=True,short_region_tail_cleared=True,
                             earlier_regions_preserved=True,non_receivers_unchanged=True,time_history_unchanged=True)
            report['checks'].append(entry)
            print('PASS',label,'exact fields/history/checkpoints, restart in both directions, buffer reuse',flush=True)
        report['status']='passed'
    except Exception as error:
        report['status']='failed';report['error']=str(error);raise
    finally:
        path=v.HERE/'shared_buffer_verification.json'
        path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('REPORT',path,flush=True);print('BUILD',v.BUILD,flush=True)

if __name__=='__main__':
    main()

