"""Prepare one protected-copy steady side-heated coarse-Sq-only test."""
from pathlib import Path
import hashlib, json, re, difflib, shutil, math

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'运行脚本/run_cases/Multiblock_SideHeated_Ra1e8_N512_flow_inverse_T1000_20261027'
REMOTE='/data2/XLLi/Multiblock/data2D/SideHeatedCell/diagnosis_N512_flow_inverse_20261027'
MASTER=ROOT/'多块网格/2DRBOpenaccMultiblock.F90'
raw=MASTER.read_bytes(); text=raw.decode('utf-8-sig').replace('\r\n','\n')
assert not OUT.exists(), 'Never overwrite an existing test'
OUT.mkdir(parents=True)
for d in ['source','source_master','results','host_reference/MB512_flow']:(OUT/d).mkdir(parents=True)
def write(path,value):path.write_text(value,encoding='utf-8',newline='\n')
def sha(data):return hashlib.sha256(data).hexdigest()
def once(s,a,b):
    assert s.count(a)==1,a
    return s.replace(a,b)

s=once(text,'!#define steadyFlow\n','#define steadyFlow\n')
s=once(s,'#define unsteadyFlow\n','!#define unsteadyFlow\n')
for a,b in [('#define RayleighBenardCell\n','!#define RayleighBenardCell\n'),
            ('#define HorizontalWallsConstT\n','!#define HorizontalWallsConstT\n'),
            ('#define VerticalWallsAdiabatic\n','!#define VerticalWallsAdiabatic\n'),
            ('!#define SideHeatedCell\n','#define SideHeatedCell\n'),
            ('!#define HorizontalWallsAdiabatic\n','#define HorizontalWallsAdiabatic\n'),
            ('!#define VerticalWallsConstT\n','#define VerticalWallsConstT\n')]:s=once(s,a,b)
s=once(s,'nx=1024, ny=1024','nx=512, ny=512')
for side,a,b in [('Left',128,64),('Right',129,65),('Bottom',128,64),('Top',129,65)]:
    s=once(s,f'fineLayerCells{side}={a}',f'fineLayerCells{side}={b}')
s=once(s,'Rayleigh=1.0d7','Rayleigh=1.0d8')
s=once(s,'Prandtl=0.7d0','Prandtl=0.71d0')
s=once(s,'itc_max=ceiling(20000000.0d0/dble(refineRatio))*refineRatio',
    'itc_max=ceiling(1000.0d0*timeUnit/dble(refineRatio))*refineRatio')
s=once(s,'outputSnapshotInterval=10.0d0','outputSnapshotInterval=0.5d0')
control=s
old='SqCoarse=1.0d0/(0.5d0+(1.0d0/Sq-0.5d0)/dxCoarse)'
new='SqCoarse=8.0d0*(2.0d0-SnuCoarse)/(8.0d0-SnuCoarse)'
s=once(s,old,new)
label="        write(k,*) 'Fine tauf, Snu, Sq, Qk, Qnu, paraA:',tauf,Snu,Sq,Qk,Qnu,paraA"
s=once(s,label,label+"\n        write(k,*) '[FLOW_ONLY] coarse Snu,Sq,Qk,Qnu,Lambda:', &\n            SnuCoarse,SqCoarse,QkCoarse,QnuCoarse,(1d0/SnuCoarse-0.5d0)*(1d0/SqCoarse-0.5d0)")
assert 'QeCoarse' not in s and 'loadInitField=0' in s
assert s.count('QnuCoarse=1.0d0/(0.5d0+(1.0d0/Qnu-0.5d0)/dxCoarse)')==1
lib,count=re.subn(r'^    program main\b.*?^    end program main\b','',s,flags=re.M|re.S)
assert count==1
write(OUT/'source/solver.F90',s);write(OUT/'source/library.F90',lib)
write(OUT/'source/original_scaling_control_not_submitted.F90',control)
(OUT/'source_master/2DRBOpenaccMultiblock.F90').write_bytes(raw)
write(OUT/'source/master_to_case.diff',''.join(difflib.unified_diff(text.splitlines(True),s.splitlines(True),fromfile='latest_protected_master',tofile='N512_flow_inverse')))
write(OUT/'source/numerical_change.diff',''.join(difflib.unified_diff(control.splitlines(True),s.splitlines(True),fromfile='same_case_original_scaling',tofile='coarse_Sq_inverse_only')))
HELP=ROOT/'运行脚本/run_cases/Multiblock_coarse_magic_20261021'
shutil.copyfile(HELP/'verify.F90',OUT/'verify.F90')
shutil.copyfile(HELP/'check_gpu_smoke.py',OUT/'check_gpu_smoke.py')
tau=.5+.1*512*math.sqrt(3*.71/1e8)
sn=1/tau;sc=1/(.5+(tau-.5)/2)
sq=8*(2-sn)/(8-sn);qc=3-math.sqrt(3);qn=4*math.sqrt(3)-6
down=lambda q:1/(.5+(1/q-.5)/2)
fine=[tau,sn,sq,qc,qn,20*math.sqrt(3)*(tau-.5)/(3*.71)-4]
coarse=[sc,8*(2-sc)/(8-sc),down(qc),down(qn)]
assert abs((1/sc-.5)*(1/coarse[1]-.5)-3/16)<1e-14
cfg=dict(remote_dir=REMOTE,name='MB512_flow',node='node05',Ra=1e8,Pr=.71,Ma=.1,
    N_fine_equivalent=512,N_coarse_equivalent=256,ratio=2,skin=2,layers=[64,65,64,65],
    center_intervals=[192,192],coarse_array=[197,197],loadInitField=0,
    target_tff_cap=1000,epsU=1e-7,epsT=1e-7,steady_check_steps=2000,
    history_snapshot_interval_tff=.5,checkpoint_tecplot_interval_tff=100,
    method='Only coarse flow Sq inverted from magic Lambda=3/16; all temperature parameters use original scaling',
    physical_initial_condition='zero velocity, linear side-heated temperature, no perturbation',
    expected_fine=fine,expected_coarse=coarse,source_sha256=sha(s.encode()),library_sha256=sha(lib.encode()),
    master_sha256=sha(raw),reference=dict(NuVolAvg=30.2225,Nu_hot=30.2423,Nu_middle=30.2223,ReVolRMS=537.342),
    walltime='12:00:00',automation=False)
write(OUT/'campaign.json',json.dumps(cfg,indent=2)+'\n')
protected=['多块网格/2DRBOpenaccMultiblock.F90','均匀网格/2DRBOpenmp.F90','均匀网格/3DRBOpenmp.F90','均匀网格/3DRBOpenacc.F90','均匀网格/3DRBOpenaccMpi.F90']
write(OUT/'protected_sources.json',json.dumps({r:sha((ROOT/r).read_bytes()) for r in protected},indent=2,ensure_ascii=False)+'\n')
write(OUT/'README.md','''# 稳态Ra=1e8，细512／粗256，仅流场反算

2026-10-27用户指定的新算例。由提交时最新多块主程序生成独立副本，主程序和既有结果不改。
侧壁加热、Ra=1e8、Pr=0.71、Ma=0.1；从零速度与原线性温度初场启动，无额外扰动、不加载检查点。
细等效512、粗等效256、r=2、细层64/65/64/65；中心192×192粗间隔、含接收层197×197数组。
只将粗格流场Sq改为由SnuCoarse通过Lambda=3/16反算；se=sepsilon=Snu绑定保持，剪切率原始缩放，温度Qk/Qnu/paraA和平衡态均不改。
温度及流场非平衡矩接口继续使用对应区域实际松弛率；碰撞、推进、边界、插值均保持原实现。不加PR/RR/H-RR。
epsU=epsT=1e-7，每2000细步检查；最多1000 TFF，达到阈值提前结束。Nu/Re与快照0.5 TFF，检查点和Tecplot100 TFF，正常结束补终点。
只提交该算例，不重复原始缩放或联合反算对照。稳态量以FinalNuRe终点为准；退出0、达到时限、无非数和稳态收敛分别判断。
当前AGENTS要求检查两节点。选择实时GPU负载较低的node05；保存检查记录，不取消或改动已有任务；轻量既有GPU进程保留，启动时核查空闲显存。
PBS先执行实际512网格40步CPU/P100对照，再启动正式求解器；节点未提供python3时，轻量校验通过SSH在master运行，不在master运行求解器。
''')
assert MASTER.read_bytes()==raw
print(json.dumps(cfg,indent=2))
