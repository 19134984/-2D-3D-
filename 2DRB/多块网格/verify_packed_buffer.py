"""Compare interface-packed histories with the previous full-coarse-buffer solver.

Usage: python verify_packed_buffer.py PATH_TO_PRE_CHANGE_SOURCE
Temporary bounds-checked gfortran OpenACC host runs; not GPU validation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
import verify_direct_exchange as d
import verify_multiblock as v
import verify_simplified as s

REGIONS=d.REGIONS

def driver(steps):
    body=d.driver(True,steps).replace('    write(77) coarseBuffer\n','')
    body=body.replace('    close(77)',"""    close(77)
    open(unit=78,file='history.bin',form='unformatted',access='stream',status='replace')
    write(78) coarseBuffer
    close(78)
    open(unit=78,file='layout.dat',status='replace')
    write(78,*) dxCoarse,interfaceSkin,centerIntegrationBox
""",1)
    for i,region in enumerate(REGIONS):
        cap=region.capitalize()
        line=f'    write(78,*) size(f_{region},1),size(f_{region},2),xOffset{cap},yOffset{cap}\n'
        if i: line='    if(refineRatio>1) then\n'+line+'    endif\n'
        # Insert immediately after prior generated layout lines.
        marker='    call exit_data_2d_openacc()'
        body=body.replace(marker,line+marker,1)
    body=body.replace('    call exit_data_2d_openacc()','    close(78)\n    call exit_data_2d_openacc()',1)
    return body

def expected_packets(folder):
    lines=(folder/'layout.dat').read_text().splitlines()
    dx,skin,xmin,xmax,ymin,ymax=map(float,lines[0].split())
    blocks=[list(map(float,line.split())) for line in lines[1:]]
    if len(blocks)==1: return np.zeros((1,1))
    ni,nj,x0,y0=blocks[0]; ni,nj=int(ni),int(nj)
    history=np.fromfile(folder/'history.bin',dtype='<f8').reshape((ni,nj,20,3),order='F')
    result=[]
    for wi,hj,x,y in blocks[1:]:
        for j in range(1,int(hj)+1):
            for i in range(1,int(wi)+1):
                px,py=x+i-.5,y+j-.5
                if not(xmin<px<xmax and ymin<py<ymax): continue
                qx,qy=(px-x0)/dx+.5,(py-y0)/dx+.5
                if abs(qx-round(qx))<1e-12 and abs(qy-round(qy))<1e-12:
                    packet=history[round(qx)-1,round(qy)-1,:,:].copy()
                else:
                    sx=max(1+int(skin),min(int(np.floor(qx))-1,ni-int(skin)-3))
                    sy=max(1+int(skin),min(int(np.floor(qy))-1,nj-int(skin)-3))
                    def weights(q):
                        return [np.prod([(q-b)/(a-b) for b in range(4) if b!=a]) for a in range(4)]
                    wx,wy=weights(qx-sx),weights(qy-sy)
                    packet=np.zeros((20,3))
                    for iy in range(4):
                        for ix in range(4):
                            packet+=wx[ix]*wy[iy]*history[sx+ix-1,sy+iy-1,:,:]
                result.append(packet)
    return np.concatenate(result,axis=0)

def packets(folder):
    ntime=3 if len((folder/'layout.dat').read_text().splitlines())>2 else 1
    return np.fromfile(folder/'history.bin',dtype='<f8').reshape((-1,ntime),order='F')

def polynomial_driver():
    body=d.polynomial_driver()
    start=body.index('    coarseBuffer=0.0d0')
    end=body.index('    rho_left=-123.0d0',start)
    body=body[:start]+"""    f_coarse=0.0d0
    g_coarse=0.0d0
    u_coarse=0.0d0
    v_coarse=0.0d0
    Fx_coarse=0.0d0
    Fy_coarse=0.0d0
    !$acc update device(f_coarse,g_coarse,u_coarse,v_coarse,Fx_coarse,Fy_coarse) async(1)
    do k=0,2
        tau=dble(k-1)
        do j=1,size(f_coarse,2)
            y=(yOffsetCoarse+(j-0.5d0)*dxCoarse)/dble(ny)
            do i=1,size(f_coarse,1)
                x=(xOffsetCoarse+(i-0.5d0)*dxCoarse)/dble(nx)
                rho_coarse(i,j)=1.0d0+0.01d0*x**3*y**3*(1.0d0+tau+tau**2)
                T_coarse(i,j)=0.02d0*x**3+0.03d0*y**3+0.04d0*tau**2
            enddo
        enddo
        !$acc update device(rho_coarse,T_coarse) async(1)
        call save_coarse_buffer(k)
        !$acc wait(1)
    enddo
"""+body[end:]
    return body

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('baseline',type=Path)
    args=ap.parse_args()
    old=args.baseline.read_text(encoding='utf-8-sig')
    new=v.SOURCE.read_text(encoding='utf-8-sig')
    assert 'coarseBuffer(:,:,:,:)' in old and 'coarseBuffer(:,:)' in new
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
            for tag,source in [('before',old),('after',new)]:
                folder,exe=v.compile_source(label+'_'+tag,s.short_settings(source,**settings),driver(40),**dims)
                v.run([exe],folder); folders.append(folder)
            before,after=folders
            err=s.compare_values(d.canonical(before),d.canonical(after),label+' fields')
            packeterr=s.compare_values(expected_packets(before),packets(after),label+' packets')
            names=['NuRe_2DOpenaccMultiblock.dat']
            if steady: names+=['Convergence_2DOpenaccMultiblock.dat']
            diagnostics={name:s.compare_values(np.loadtxt(before/name),np.loadtxt(after/name),label+' '+name) for name in names}
            snapshoterr=s.compare_snapshots(next(before.glob('*Snapshot-*.bin')),next(after.glob('*Snapshot-*.bin')))
            restored,exe=v.compile_source(label+'_read_old',s.short_settings(new,restart=True,**settings),driver(0),**dims)
            d.copy_checkpoint(before,restored); v.run([exe],restored)
            s.compare_values(d.canonical(before),d.canonical(restored),label+' restored fields',exact=True)
            conversionerr=s.compare_values(expected_packets(before),packets(restored),label+' converted history')
            half=(20//ratio)*ratio
            resume_errors={}
            for tag,source in [('old',old),('new',new)]:
                split,exe=v.compile_source(label+'_'+tag+'_half',s.short_settings(source,**settings),driver(half),**dims)
                v.run([exe],split)
                _,exe=v.compile_source(label+'_'+tag+'_resume',s.short_settings(new,restart=True,**settings),driver(40-half),**dims)
                v.run([exe],split)
                resume_errors[tag]=s.compare_values(d.canonical(after),d.canonical(split),label+' '+tag+' resume',exact=(tag=='new'))
                s.compare_values(packets(after),packets(split),label+' '+tag+' resume history',exact=(tag=='new'))
            constant,exe=v.compile_source(label+'_constant',s.short_settings(new,**settings),s.constant_driver(),**dims)
            v.run([exe],constant)
            actual,expected=np.loadtxt(constant/'constant_diagnostics.dat')
            s.compare_values(expected,actual,label+' constant diagnostics')
            entry=dict(case=label,fields_max_abs_error=err,packets_max_abs_error=packeterr,
                       diagnostics=diagnostics,snapshot_max_abs_error=snapshoterr,
                       geometry_and_integration_areas_exact=True,old_checkpoint_fields_exact=True,
                       old_history_conversion_max_abs_error=conversionerr,restart_max_abs_errors=resume_errors,
                       old_history_bytes=(before/'history.bin').stat().st_size,
                       packed_history_bytes=(after/'history.bin').stat().st_size,
                       old_checkpoint_bytes=d.c.checkpoint(before).stat().st_size,
                       new_checkpoint_bytes=d.c.checkpoint(after).stat().st_size)
            if ratio>1:
                poly,exe=v.compile_source(label+'_polynomial',s.short_settings(new,**settings),polynomial_driver(),**dims)
                v.run([exe],poly)
                count,error=np.loadtxt(poly/'polynomial.dat')
                entry.update(polynomial_receivers=int(count),polynomial_max_abs_error=float(error),non_receiver_nodes_unchanged=True)
            report['checks'].append(entry)
            print('PASS',label,'fields, all packet components, old/new restart, area, polynomial',flush=True)
        report['analytical_statistics']=s.verify_statistics(new)
        report['status']='passed'
    except Exception as error:
        report['status']='failed'; report['error']=str(error); raise
    finally:
        path=v.HERE/'packed_buffer_verification.json'
        path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('REPORT',path,flush=True); print('BUILD',v.BUILD,flush=True)

if __name__=='__main__':
    main()

