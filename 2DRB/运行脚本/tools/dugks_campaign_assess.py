"""Independent DUGKS host/GPU, finite-checkpoint and statistical assessments."""
from array import array
from pathlib import Path
import argparse
import json
import math
import re
import struct
import sys


def record(stream):
    prefix=stream.read(4)
    if not prefix:
        return None
    if len(prefix)!=4:
        raise ValueError("Truncated Fortran record marker")
    length=struct.unpack("<i",prefix)[0]
    if length<0:
        raise ValueError("Unsupported split Fortran record")
    data=stream.read(length)
    if len(data)!=length or stream.read(4)!=prefix:
        raise ValueError("Truncated Fortran record")
    return data


def doubles(data):
    result=array("d")
    result.frombytes(data)
    if sys.byteorder!="little":
        result.byteswap()
    return result


def smoke(folder):
    cfg=json.loads((folder/"case.json").read_text())
    result=dict(case=cfg["name"],passed=True,source_sha256=cfg["source_sha256"],
                smoke_sha256=cfg["smoke_sha256"],total_steps=45,measured_steps=40,fields={})
    with (folder/"smoke/host/smoke.bin").open("rb") as host,(folder/"smoke/gpu/smoke.bin").open("rb") as gpu:
        hh,gh=record(host),record(gpu)
        assert hh==gh and struct.unpack("<2i",hh)==(cfg["nx"],cfg["ny"])
        for name in ["f_tilde","g_tilde","u","v","T","rho"]:
            a,b=doubles(record(host)),doubles(record(gpu))
            assert len(a)==len(b)
            largest=0.0
            for x,y in zip(a,b):
                if not math.isfinite(x) or not math.isfinite(y):
                    raise ValueError("Nonfinite smoke field "+name)
                largest=max(largest,abs(x-y))
            result["fields"][name]=dict(count=len(a),max_absolute_error=largest,passed=largest<=1e-10)
            result["passed"]=result["passed"] and largest<=1e-10
            if name=="rho":
                result["gpu_relative_mass_error"]=abs(sum(b)/(cfg["nx"]*cfg["ny"])-1)
                if cfg.get("restart"):
                    # 续算场已包含种子累计舍入误差；另检验本次短程新增质量漂移。
                    seed=json.loads((folder/"results/restart_preparation.json").read_text())
                    initial=seed["mean_density"]
                    result["mass_reference"]="verified checkpoint mean density"
                    result["initial_relative_mass_error"]=abs(initial-1)
                    result["gpu_incremental_mass_error"]=abs(sum(b)/(cfg["nx"]*cfg["ny"])-initial)/abs(initial)
                    assert result["gpu_incremental_mass_error"]<1e-10
                else:
                    assert result["gpu_relative_mass_error"]<1e-10
        assert record(host) is None and record(gpu) is None
    for device in ["host","gpu"]:
        log=(folder/("smoke/"+device+"/solver.stdout")).read_text(errors="replace")
        elapsed=float(re.search(r"SMOKE_ELAPSED_SECONDS\s+(\S+)",log).group(1))
        result[device+"_seconds_40_steps"]=elapsed
        result[device+"_MLUPS"]=cfg["nx"]*cfg["ny"]*40/max(elapsed,1e-12)/1e6
    result["estimated_1000_tff_hours"]=cfg["step_cap"]*result["gpu_seconds_40_steps"]/40/3600
    result["timing_limit"]="Short warm-start shared-GPU estimate; excludes I/O and long-run contention."
    return result


def integrate(rows,start,end,tolerance):
    if not rows or rows[0][1]>start+tolerance or rows[-1][1]<end-tolerance:
        return dict(complete=False,start=start,end=end)
    nu=rr=duration=0.0
    for left,right in zip(rows,rows[1:]):
        t0,t1=left[1],right[1]
        a,b=max(start,t0),min(end,t1)
        if b<=a:
            continue
        frac0,frac1=(a-t0)/(t1-t0),(b-t0)/(t1-t0)
        n0=left[2]+(right[2]-left[2])*frac0
        n1=left[2]+(right[2]-left[2])*frac1
        r0=left[3]+(right[3]-left[3])*frac0
        r1=left[3]+(right[3]-left[3])*frac1
        nu+=(b-a)*(n0+n1)/2
        rr+=(b-a)*(r0*r0+r0*r1+r1*r1)/3
        duration+=b-a
    if duration<=0:
        return dict(complete=False,start=start,end=end)
    return dict(complete=abs(duration-(end-start))<=2*tolerance,start=start,end=end,
                covered_duration=duration,Nu=nu/duration,Re=math.sqrt(rr/duration))


def finished(folder):
    cfg=json.loads((folder/"case.json").read_text())
    result=folder/"results"
    status=int((result/"run.status").read_text())
    logs="\n".join((result/name).read_text(errors="replace") for name in ["solver.stdout","solver.stderr"])
    report=dict(case=cfg["name"],source_sha256=cfg["source_sha256"],solver_exit_status=status,
        diagnostic_error=bool(re.search(r"(?i)\b(nan|infinity|cuda_error|error\s+stop|error:)\b",logs)))
    lines=(result/"reloadFile2DOpenaccDUGKS-latest.meta").read_text().splitlines()
    magic="DUGKS2D_D2Q9_D2Q9_V1"
    assert lines[0].strip()==magic
    final_itc,final_tf=map(float,lines[2].split())
    report.update(final_step=int(final_itc),final_tff=final_tf,reached_horizon=final_itc>=cfg["step_cap"])
    with (result/lines[1].strip()).open("rb") as stream:
        header=stream.read(32+5*4+20*8)
        assert header[:32].rstrip()==magic.encode()
        nx,ny,casecode,flowcode,periodic=struct.unpack("<5i",header[32:52])
        assert (nx,ny)==(cfg["nx"],cfg["ny"]) and periodic==0
        assert (casecode,flowcode)==((1,1) if cfg["case_type"]=="side" else (2,2))
        state=stream.read(7*4+2*8)
        nstep,nrow,nsnap,nplt,nreload,nconv,last,eu,et=struct.unpack("<7i2d",state)
        assert nstep==int(final_itc)
        count=nx*ny
        mass=None
        for name in ["f","g"]:
            values=doubles(stream.read(9*count*8))
            assert len(values)==9*count and all(math.isfinite(x) for x in values)
            if name=="f":
                mass=array("d",[0.0])*count
                for k,x in enumerate(values):
                    mass[k%count]+=x
        assert min(mass)>0
        references=doubles(stream.read(3*count*8))
        assert len(references)==3*count and all(math.isfinite(x) for x in references)
        cached=doubles(stream.read(10*nrow*8))
        assert len(cached)==10*nrow and all(math.isfinite(x) for x in cached)
        assert stream.read(32).rstrip()==magic.encode() and stream.read(1)==b""
    report.update(checkpoint_finite=True,minimum_density=min(mass),
                  relative_mass_error=abs(sum(mass)/count-1),history_rows=nrow)
    rows=[]
    for line in (result/"NuRe_2DOpenaccDUGKS.dat").read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        row=list(map(float,line.split()))
        assert len(row)==10 and all(math.isfinite(x) for x in row)
        assert not rows or row[1]>rows[-1][1]
        rows.append(row)
    assert len(rows)==nrow
    report["last_sample"]=dict(time_tff=rows[-1][1],Nu=rows[-1][2],Re=rows[-1][3])
    final_lines=(result/"FinalDiagnostics_2DOpenaccDUGKS.dat").read_text().splitlines()
    index=next(i for i,line in enumerate(final_lines) if line.strip()=="NuVolAvg ReVolAvg Nu_global Nu_hot Nu_cold Nu_middle")
    values=list(map(float,final_lines[index+1].split()))
    assert len(values)==6 and all(math.isfinite(x) for x in values)
    report["final_state_diagnostics"]=dict(zip(
        ["NuVolAvg","ReVolAvg","Nu_global","Nu_hot","Nu_cold","Nu_middle"],values))
    report["integration_succeeded"]=status==0 and not report["diagnostic_error"]
    if cfg["case_type"]=="side":
        report.update(errorU=eu,errorT=et,consecutive_passes=nconv,
                      steady_converged=nconv>=3 and eu<1e-7 and et<1e-7)
    else:
        tolerance=2*cfg["dt"]/(cfg["nx"]/cfg["Uff"])
        windows=[integrate(rows,a,b,tolerance) for a,b in cfg["average_windows"]]
        report["windows"]=windows
        if all(w["complete"] for w in windows):
            diff={name:100*abs(windows[0][name]-windows[1][name])/abs(windows[1][name])
                  for name in ["Nu","Re"]}
            report["half_window_difference_percent"]=diff
            report["half_window_check_passed"]=all(x<1 for x in diff.values())
            report["reference_error_percent"]={
                name:100*(windows[1][name]/cfg["reference_"+name]-1) for name in ["Nu","Re"]}
    report["passed"]=report["integration_succeeded"] and (
        report.get("steady_converged",False) if cfg["case_type"]=="side" else report["reached_horizon"])
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--case",type=Path,required=True)
    args=parser.parse_args()
    try:
        result=finished(args.case)
    except Exception as exc:
        result=dict(passed=False,error=str(exc))
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result["passed"] else 87)
