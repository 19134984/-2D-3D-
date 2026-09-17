"""Regression: unchanged host trajectory and GPU/host field agreement after 40 fine steps."""
from pathlib import Path
import sys, json, hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'多块网格'))
import verify_multiblock as v
OUT=ROOT/'运行脚本/Multiblock_gpu_debug'
DRIVER='''program main
    use commondata
    use openacc
    implicit none
    integer :: step
    call acc_init(acc_device_default)
    call initial()
    call enter_data_2d_openacc()
    do step=1,20
        call advance_multiblock()
    enddo
    call update_host_all(.true.)
    call calNuRe()
    if (.not.all(ieee_is_finite(fStorage))) error stop 'Nonfinite flow population'
    if (.not.all(ieee_is_finite(gStorage))) error stop 'Nonfinite thermal population'
    open(unit=77,file='state.bin',form='unformatted',access='stream',status='replace')
    write(77) fStorage,gStorage,uStorage,vStorage,TStorage,rhoStorage,pStorage
    close(77)
    call exit_data_2d_openacc()
    print *, 'REGRESSION_FINE_STEPS',itc
end program main'''

if __name__=='__main__':
    old=(ROOT/'运行脚本/Multiblock_sideheated_20260916_v2/Ra1e6/ratio2/sideheated_N256_20260916_v2/source/solver.F90').read_text(encoding='utf-8')
    fixed=(OUT/'fixed_case.F90').read_text(encoding='utf-8')
    states=[]
    for name,s in [('before',old),('after',fixed)]:
        folder,exe=v.compile_source('fix_'+name,s,DRIVER,n=256,ra=1e6)
        print(name,v.run([exe],folder,timeout=120).strip(),flush=True)
        states.append(np.fromfile(folder/'state.bin',dtype='<f8'))
        if name=='after':
            (OUT/'host_state.bin').write_bytes((folder/'state.bin').read_bytes())
            (OUT/'gpu_regression.F90').write_text((folder/'solver.F90').read_text(encoding='utf-8'),encoding='utf-8',newline='\n')
    assert all(np.all(np.isfinite(s)) for s in states)
    assert np.array_equal(*states),np.max(np.abs(states[0]-states[1]))
    report=dict(fine_steps=40,n=256,ra=1e6,host_before_after_bitwise_equal=True,
                state_values=len(states[0]),source_sha256=hashlib.sha256(v.SOURCE.read_bytes()).hexdigest())
    (OUT/'host_regression.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
