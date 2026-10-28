"""Stdlib-only CPU/GPU field comparison for the complete coupled solver."""
from pathlib import Path
import json
import math
import struct
import sys

root=Path(__file__).resolve().parent


def read_snapshot(path):
    fields=[]
    with path.open('rb') as stream:
        for _ in range(4):
            size=struct.unpack('<i',stream.read(4))[0]
            values=struct.unpack('<'+'d'*(size//8),stream.read(size))
            assert size==struct.unpack('<i',stream.read(4))[0]
            assert all(math.isfinite(x) for x in values)
            fields.append(values)
        assert stream.read()==b''
    return fields


reports={}
for strategy in [1,2]:
    case=root/('GPU_P'+str(strategy))
    host=list((case/'host_reference').glob('*Snapshot*.bin'))
    gpu=list((case/'results').glob('*Snapshot*.bin'))
    assert len(host)==len(gpu)==1,(host,gpu)
    a=read_snapshot(host[0]); b=read_snapshot(gpu[0])
    errors=[max(abs(x-y) for x,y in zip(p,q)) for p,q in zip(a,b)]
    assert max(errors)<1e-8,errors
    reports['P'+str(strategy)]={'max_abs_u_v_T_rho':errors,'pass':True}
out=root/'GPU_validation/results/gpu_comparison.json'
out.write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports,indent=2))
