"""汇集冻结结果；只生成报告，不运行或修改求解器。"""
from pathlib import Path
import json,csv,math
P=Path(__file__).resolve().parent
REPO=P.parent.parent
C=REPO/'运行脚本/run_cases'
stats=json.loads((P/'independent_statistics.json').read_text())
rows=[]
REF={6:(6.30,279),7:(11.37,968),8:(25.25,3662),9:(53.51,15101),10:(93.88,61797)}

def add(group,exp,N,method,status,tff=None,nu=None,re=None,window='',note='',evidence='',uniform=False,layers=None):
    if layers is None:layers=[N//8,N//8+1]*2
    r=dict(id='T{:02d}'.format(len(rows)+1),group=group,geometry='稳态侧壁加热' if group=='稳态基准' else '非稳态RB',
           Ra=10**exp,N_fine_equivalent=N,grid='单块均匀' if uniform else '五块粗细网格',ratio=1 if uniform else 2,
           uniform_coarse_equivalent=None if uniform else N//2,
           fine_layers=None if uniform else layers,
           center_coarse_intervals=None if uniform else (N-layers[0]-layers[1]+1)//2,
           algorithm=method,status=status,tff=tff,window=window,Nu=nu,Re=re,
           Nu_error_percent=None,Re_error_percent=None,note=note,evidence=evidence)
    if nu is not None and group!='稳态基准':r['Nu_error_percent']=100*abs(nu/REF[exp][0]-1)
    if re is not None and group!='稳态基准':r['Re_error_percent']=100*abs(re/REF[exp][1]-1)
    rows.append(r)
    return r

orig='D2Q9流场TRT＋D2Q5 Luo温度；粗区各松弛时间平移量原始缩放'
magic='流场TRT：黏性缩放＋魔法参数反算；温度仍原始缩放'
thermal='流场仍原始缩放；温度扩散缩放＋Luo四阶关系反算（粗区MRT）'
joint='联合反算：流场TRT魔法参数＋温度Luo四阶关系（粗区MRT）'
reg=lambda name:name+'流场正则化；温度D2Q5 Luo原始缩放'
useg='D2Q9流场TRT魔法参数反算＋D2Q5 EnableUseG温度'
cde='双D2Q9、双BGK的LBM-CDE反向设计'

side=json.loads((C/'Multiblock_data2D_SideHeatedCell_N1024_R2_20261019/results_audit_7265.json').read_text())
for e,x in zip((6,7,8),side['cases']):
    r=add('稳态基准',e,1024,orig,'收敛（epsU、epsT均1e-7）',x['t_ff'],x['final']['NuVolAvg'],x['final']['ReVolRMS'],
          '最终收敛时刻',evidence='Multiblock_data2D_SideHeatedCell_N1024_R2_20261019/results_audit_7265.json')
    r['Nu_error_percent']=x['relative_error_percent']['NuVolAvg'];r['Re_error_percent']=x['RMS_vs_table_percent']
    r['four_metric_values']={k:x['final'][k] for k in ['NuVolAvg','Nu_hot','Nu_middle','ReVolRMS']}
    r['four_metric_errors']={k:x['relative_error_percent'][k] for k in ['NuVolAvg','Nu_hot','Nu_middle']}
    r['four_metric_errors']['ReVolRMS']=x['RMS_vs_table_percent']
    r['note']='Pr=0.71；对照论文表2的1025²行；Re数值采用RMS比较，论文式(34)定义疑点仍保留。'

ev='Multiblock_data2D_RB_Pr0p7_T1000_R2_20261019/postprocess_20261020_0905/结果汇总.md'
add('原始RB基准',6,128,orig,'完成，但选中上下双涡分支',1000.340531,4.42324812566,218.389036735,'500–1000',note='被加扰动重算替代；不能作为最终Ra6基准。',evidence=ev)
add('原始RB基准',7,256,orig,'完成1000 TFF',1000.439763,11.368352446,965.386059598,'500–1000',evidence=ev)
for e,N,t in [(8,512,18.508135622),(9,1024,14.503106420),(10,2048,14.001420610)]:
    add('原始RB基准',e,N,orig,'无效场值退出',t,evidence=ev)
ev='Multiblock_data2D_RB_Pr0p7_T1000_R2_20261019/diagnosis_20261020/诊断结论与松弛参数核对.md'
add('原始RB基准',6,128,orig+'；初始温度加入1e-3正弦扰动','完成200 TFF诊断',200.2683746251514,6.30872037214,280.413654429,'100–200',evidence=ev)
add('原始RB基准',6,128,orig+'；初始温度加入1e-3正弦扰动','从零重算完成1000 TFF',1000.340531,6.315109604,280.003470570,'500–1000',evidence='Multiblock_RB_Ra1e6_seed_T1000_20261020/postprocess_20261020_Ra6rerun_Ra7/结果汇总.md')

for method,status,t,N,u,note in [
    (orig+'；从零重复','无效场值退出',18.508135622,512,False,'与原历史逐字节一致'),
    (orig+'；同一有效检查点续算','无效场值退出',18.508135622,512,False,'用于验证可重复性'),
    (orig+'；粗到细空间插值改双线性','无效场值退出',18.508135622,512,False,''),
    (orig+'；二次时间插值改线性','无效场值退出',18.508135622,512,False,''),
    ('D2Q9流场TRT＋D2Q5 Luo温度，匹配细级参数；单块无跨级缩放','完成30 TFF',30.006427076,512,True,'无粗细接口；只做短测'),
    (orig+'；加入Ra6同款初始扰动','无效场值退出',10.004397633,512,False,''),
    (magic+'；旧检查点暖启动','无效场值退出',18.007915740,512,False,'旧接口时间历史未重建；不能作为单独否定证据'),
    (magic+'；从零启动','无效场值退出',18.007915740,512,False,''),
    ('均匀网格显式匹配多块粗区原缩放的流场及温度参数','完成30 TFF',30.013192900,256,True,''),
    ('均匀网格匹配粗区参数，但流场魔法参数恢复3/16','完成30 TFF',30.013192900,256,True,'')]:
    add('缩放与插值隔离',8,N,method,status,t,note=note,evidence=ev,uniform=u)

add('网格与分区试验',8,512,orig+'；扩大四周细区','无效场值退出',16.507256094947547, layers=[128,129,128,129],evidence='Multiblock_Ra1e8_wide_fine_T30_20261020/短试验结果.md')
add('网格与分区试验',8,1024,orig+'；整体加密一倍','目标终点检查失败，不能算完成',30.006427,evidence='Multiblock_Ra1e8_refined_N1024_T30_20261020/短试验结果.md')
ev='Multiblock_Ra8_Ra9_MRTclosure_diag_20261020/诊断结论.md'
for e,N,method,t in [(8,512,thermal,18.508136),(8,512,joint,18.007916),(9,1024,magic,14.002999),(9,1024,thermal,14.503106),(9,1024,joint,14.002999)]:
    add('Luo反算闭合',e,N,method,'无效场值退出；未完成30 TFF',t,evidence=ev)

ev='Multiblock_interface_audit_20261020/检查结论.md'
for name,status in [('原对照','失败'),('接口仅抑制m1、m2能量类非平衡矩','通过该局部窗口'),('接口仅抑制m4、m6奇矩非平衡部分','失败'),('接口同时抑制能量类和奇矩非平衡部分','通过该局部窗口')]:
    add('接口模态局部诊断',8,512,orig+'；'+name,status,18.213596773,
        note='从累计161000步有效检查点推进到161520步；不是从零通过30 TFF。',evidence=ev)
add('接口与滤波短测',8,512,orig+'；只在接口抑制能量类非平衡矩','无效场值退出',19.008356,evidence=ev)
add('接口与滤波短测',8,512,orig+'；粗区m1、m2松弛率改1.2','无效场值退出',21.009235,note='体积压缩响应也改变，不能称只调无物理作用的误差项。',evidence=ev)
ev='Multiblock_interface_literature_20261020/文献方案试验结果.md'
for name,t in [('细到粗九点加权非平衡滤波',18.508136),('细到粗九点等权非平衡滤波',18.508136),('九点加权滤波＋粗区能量矩率1.2',20.509015)]:
    add('接口与滤波短测',8,512,orig+'；'+name,'无效场值退出',t,evidence=ev)

ev='Multiblock_regularization_20261020/短测结果.md'
add('正则化短测',8,512,orig+'；只刷新当前力速度','无效场值退出',18.508136,evidence=ev)
for name,status,t in [('PR2','完成30 TFF',30.013193),('RR3','完成30 TFF',30.013193),('H-RR3 σ=0.98','完成30 TFF',30.013193),('H-RR3 σ=0.90','无效场值退出',3.001319),('H-RR3 σ=0.90，块首末两层回退RR','完成30 TFF',30.013193)]:
    add('正则化短测',8,512,reg(name),status,t,note='H-RR为接口梯度简化的工程适配；短测不是文献误差。',evidence=ev)

names={'pr2':'PR2','rr3':'RR3','hrr98':'H-RR3 σ=0.98','hrr90_boundary_rr':'H-RR3 σ=0.90，块边缘RR'}
for path,x in stats.items():
    parts=Path(path).parts
    method=names.get(parts[-2],'RR3')
    if 'diagnosis_regularization_T1000_' in path:
        e=8;group='正则化1000TFF';status='完成1000 TFF'
    elif 'diagnosis_regularization_Ra9_Ra10_' in path:
        e=10 if 'Ra1e10' in path else 9;group='高Ra正则化';status='无效场值退出'
    else:
        e=8 if 'Ra1e8' in path else 9;group='加密RR3';status='运行中，目标1000 TFF' if x['status']=='running' else '完成100 TFF短测'
    long=x.get('window')==[500,1000]
    # 运行中不重复列启动阶段窗口为最终统计。
    values=(x.get('Nu'),x.get('Re_space_time_RMS')) if x['status']!='running' else (None,None)
    r=add(group,e,x['N'],reg(method),status,x.get('invalid_detected_tff',x['last_tff']),*values,
          window='500–1000' if long else ('50–100（早期诊断）' if values[0] is not None else ''),
          note='长期统计尚未完成；截至2026-10-23 13:45约569.56 TFF。' if x['status']=='running' else '',
          evidence=str(Path('文档/多块网格测试总表_20261023/server_evidence')/path))
    r['mass_change_percent']=x['mass_change_percent'];r['max_mass_drift_percent']=x['max_mass_drift_percent']
    r['quality_diagnostics']=x

ev='Multiblock_UseG_T30_20261020/短测结果.md'
add('UseG与均匀对照',8,512,'原始TRT流场（Sq仍缩放）＋EnableUseG温度','无效场值退出',18.508136,evidence=ev)
add('UseG与均匀对照',8,512,useg,'无效场值退出；未完成30 TFF',18.007916,evidence='Multiblock_UseG_magic_T30_20261021/短测结果.md')
ev='UseG_uniform_multiblock_convergence_20261021/研究结论_粗细两种均匀对照.md'
add('UseG与均匀对照',8,512,useg+'；长算任务独立重复','无效场值退出',18.007916,evidence=ev)
add('UseG与均匀对照',8,256,'D2Q9流场TRT（魔法参数3/16）＋D2Q5 EnableUseG；匹配粗级参数','完成1500 TFF，通过预设统计筛选',1500.159425,25.11785158,3643.92268867,'500–1500',evidence=ev,uniform=True)
add('UseG与均匀对照',8,512,'D2Q9流场TRT（魔法参数3/16）＋D2Q5 EnableUseG；匹配细级参数','完成1500 TFF，通过预设统计筛选',1500.321354,25.24602382,3666.89085652,'500–1500',evidence=ev,uniform=True)

early=json.loads((C/'Multiblock_LBMCDE_Ra8_20261021/completed_window_diagnostics_20261022.json').read_text())
for x in early:
    u=x['case'].startswith('uniform');N=256 if '256' in x['case'] else 512
    chi={'cde_chi0':0,'cde_chi08':.8,'cde_chi095':.95,'cde_chi099':.99,'uniform256_chi08':.8,'uniform512_chi08':.8}[x['case']]
    add('早期双BGK留档',8,N,cde+'；全域χν=χκ={}，χb=0'.format(chi),'完成100 TFF；多块热halo有遗漏，结果作废' if not u else '完成100 TFF；早期窗口',100.021424 if u and N==512 else 100.043976,
        x['Nu'],x['Re_space_time_RMS'],'50–100（早期诊断）',
        note='多块温度D2Q9实际只交换5方向，后续已修正；该批不能用作有效算法精度依据。' if not u else '单块分支没有细细halo交换，不受该遗漏直接影响；但未完成长期统计。',
        evidence='Multiblock_LBMCDE_Ra8_20261021/完整短测结果_20261022.md',uniform=u)

ev='Multiblock_LBMCDE_LevelChi_Ra8_Ra9_20261022/测试结果_20261022.md'
add('修正halo后分区χ',8,512,cde+'；细χ=0、粗χ=0.5','完成100 TFF短测',100.043976,24.21245944,4739.64562467,'50–100（早期诊断）',note='三类粗区χ均0.5，包含体黏性匹配；各分布粗细基础τ分别相同。',evidence=ev)
add('修正halo后分区χ',9,1024,cde+'；细χ=0、粗χ=0.5','无效场值退出',32.006856,evidence=ev)
ev='Multiblock_LBMCDE_Ra9_CoarseChiScan_20261022/测试结果_20261022.md'
for chi,nu,re,m in [(.8,49.42668005,14942.132816,.00522774),(.9,47.88701229,13001.552469,1.3070699)]:
    r=add('修正halo后分区χ',9,1024,cde+'；细χ=0、粗χ='+str(chi),'完成100 TFF短测',100.021424,nu,re,'50–100（早期诊断）',evidence=ev)
    r['max_mass_drift_percent']=m
latest=json.loads((C/'Multiblock_LBMCDE_Ra9_Chi08_09_T1000_20261022/results_audit_20261023.json').read_text())['cases']
for chi,k in [(.8,'Ra1e9_coarse08'),(.9,'Ra1e9_coarse09')]:
    x=latest[k];good=x['exit_status']==0;w=x.get('windows',{}).get('500-1000',{})
    r=add('修正halo后分区χ',9,1024,cde+'；细χ=0、粗χ='+str(chi),
          '完成1000 TFF，但质量严重流失，未通过验收' if good else '续算后无效场值退出',x.get('failure_tff',x['last_tff']),w.get('Nu_volume'),w.get('Re_space_time_RMS'),
          '500–1000（有严重质量缺陷的诊断值）' if good else '',note='总质量下降约26.02%，不能作为可靠基准结果。' if good else '最后历史302.564806 TFF。',
          evidence='Multiblock_LBMCDE_Ra9_Chi08_09_T1000_20261022/结果复核_20261023.md')
    r['mass_change_percent']=100*x['final_mass_relative_change'];r['max_mass_drift_percent']=100*x['max_absolute_relative_mass_drift']

(P/'全部测试清单.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
fields=['id','group','geometry','Ra','grid','N_fine_equivalent','uniform_coarse_equivalent','ratio','fine_layers','center_coarse_intervals','algorithm','status','tff','window','Nu','Re','Nu_error_percent','Re_error_percent','mass_change_percent','max_mass_drift_percent','note','evidence']
with (P/'全部测试清单.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

def fmt(x,n=4):return '—' if x is None else ('{:.'+str(n)+'f}').format(x)
def grid(r):return ('均匀' if r['ratio']==1 else '多块F')+str(r['N_fine_equivalent'])+'²'

lines=['# 本轮多块网格数值测试总表（2026-10-23）','',
'覆盖本次对话的稳态验证、RB正式算例、缩放/联合反算、接口试验、正则化、EnableUseG、双D2Q9-BGK反向设计以及CPU辅助验证。截止服务器2026-10-23 13:45的状态；正在运行的算例随后可能继续推进。',
'','## 读表约定','',
'- 多块F N²表示四周细级格距等效于全域N×N均匀网格，并非所有区域都按N²计算。粗细比均为2；中心粗级格距等效全域(N/2)²。均匀N²表示单块全域N×N、粗细比1。',
'- 默认细层参数按左/右/下/上为N/8、N/8+1、N/8、N/8+1。中心粗区间隔数为3N/8，数组边长另含接收层，为3N/8+5。扩大细区算例例外，已单列。',
'- 除稳态Pr=0.71外，RB均Pr=0.7、Ma=0.1。表中失败TFF指程序检测时刻，通常每0.5 TFF检查；不代表精确的首次失稳时刻。',
'- Re一律指RMS：稳态为空间RMS，非稳态为sqrt(时间平均(瞬时空间RMS Re的平方))；Nu为体积定义的时间平均。相对误差为100×|计算−文献|/|文献|。',
'- 完成30/100 TFF只说明通过该有限时长；50–100或100–200窗口的误差仅为早期诊断偏差。代码存在已知缺陷或质量严重流失时，即使给出数值偏差，也不能作为精度通过证据。',
'- 同一方案的短测、长测在本表分别列出；同一控制组的多次复现归入对应行或附录，避免将重复作业误当新算法。',
'','## 网格明细','',
'|细级等效N²|粗级格距等效全域网格|四周细层参数（左/右/下/上）|中心粗区间隔|粗区数组（含接收层）|',
'|---|---|---|---|---|']
for n in (128,256,512,1024,2048):
    a=n//8;b=a+1;k=3*n//8
    lines.append(f'|{n}²|{n//2}²|{a}/{b}/{a}/{b}|{k}×{k}|{k+5}×{k+5}|')
lines += ['|512²（扩大细区）|256²|128/129/128/129|128×128|133×133|','',
'## 算法名称的准确含义','',
'1. **原始缩放**：细级流场TRT、温度D2Q5 Luo特殊TRT；粗级保留同一温度平衡参数a，Snu、Sq、Qk、Qnu均满足H(s粗)=H(s细)/2，H(s)=1/s−1/2。原细级满足的魔法参数和温度四阶关系一般不再在粗级成立。',
'2. **流场魔法参数反算**：黏性仍按粗细缩放，另一组TRT率由H(Snu粗)H(Sq粗)=3/16确定；流场仍是两组率的TRT。',
'3. **温度四阶反算**：Qk按扩散缩放，Qe、Qnu由Luo四阶关系分别确定；粗级Qe通常不等于Qnu，因此温度变成MRT，细级仍保留原特殊TRT。**联合反算**指2和3同时采用。',
'4. **PR2/RR3/H-RR3**：改变全流场非流体模态的投影/递归/混合重构；温度仍是原D2Q5 Luo缩放。它们不是联合反算，也不继续使用原TRT全部矩的绑定。H-RR为记录过接口梯度简化的工程适配。',
'5. **EnableUseG**：移植均匀核心已有D2Q5温度分支，含Tu时间修正及配套接口历史；在本轮没有把它叫作Luo温度MRT。',
'6. **LBM-CDE反向设计**：流场和温度都为D2Q9、BGK。早期全域χ扫描令χν=χκ=χ、χb=0；后期分区χ令细级三种χ=0、粗级三种χ同值以保持物理参数匹配。这两种参数试验并不完全相同。',
'','## 文献对照值','',
'本次重新读取仓库原始PDF表格核对。RB参照Zhang、Zhou、Sun (2017) Table 1：Ra=10⁶/10⁷/10⁸/10⁹/10¹⁰的Nu为6.30/11.37/25.25/53.51/93.88，Re为279/968/3662/15101/61797。',
'','稳态参照ISLBM文献Table 2的1025²行。文献式(34)写为空间平均速度大小，而表中Re与我们的RMS接近；本表按当前约定列RMS数值比较，保留该定义疑点，不判定作者写错。',
'','## 稳态四项指标','',
'|Ra|Nu体积（误差%）|Nu热壁（误差%）|Nu中线（误差%）|Re空间RMS（与表值差异%）|',
'|---|---|---|---|---|']
for r in rows[:3]:
    vals=r['four_metric_values'];errs=r['four_metric_errors']
    lines.append('|10^{}|{}|'.format(round(math.log10(r['Ra'])),'|'.join(fmt(vals[k],6)+'（'+fmt(errs[k],4)+'%）' for k in ['NuVolAvg','Nu_hot','Nu_middle','ReVolRMS'])))
lines += ['','## 全部耦合算例与隔离试验','']
last=None
for r in rows:
    if r['group']!=last:
        lines += ['','### '+r['group'],'', '|编号|Ra|网格|算法与改动|状态 / 完成或检测TFF|Nu / Re|Nu / Re相对偏差%|窗口|', '|---|---|---|---|---|---|---|---|']
        last=r['group']
    lines.append('|'+ '|'.join([r['id'],'10^'+str(round(math.log10(r['Ra']))),grid(r),r['algorithm'],r['status']+'；'+fmt(r['tff'],6),fmt(r['Nu'],6)+' / '+fmt(r['Re'],3),fmt(r['Nu_error_percent'],4)+' / '+fmt(r['Re_error_percent'],4),r['window'] or '不足统计窗'])+'|')
lines += ['','## 长算质量诊断','',
'|方案|Ra / 网格|质量最大相对漂移%|Nu半窗差%|Re半窗差%|冷热壁与体积Nu最大差%|','|---|---|---|---|---|---|']
for r in rows:
    if r['group']=='正则化1000TFF':
        x=r['quality_diagnostics'];lines.append('|{}|10^8 / {}|{}|{}|{}|{}|'.format(r['algorithm'].split('；')[0],grid(r),fmt(x['max_mass_drift_percent'],6),fmt(x['Nu_half_difference_percent']),fmt(x['Re_half_difference_percent']),fmt(x['wall_volume_difference_percent'])))
lines += ['','四组正则化长算均完成1000 TFF，Nu/Re偏差均小于1%；PR2及H-RR σ=0.98的Re两半窗口差仍略超过1%。短时通过、统计筛选、文献吻合及网格无关性仍是不同结论。',
'','分区χ的Ra9结果尤其要注意：χ粗=0.8完成100 TFF后，在303.064913 TFF退出；χ粗=0.9完成1000.214236 TFF，但质量净下降26.019140%，500–1000诊断Nu=45.485154、Re=9669.4831，偏离文献14.9969%、35.9679%。这组不能作为可靠基准。',
'','加密RR3的Ra9、F2048²长算仍在运行，收集时为569.557788 TFF，总质量相对初始变化约+0.00001156%。未覆盖500–1000，不给出最终误差。Ra8/F1024²和Ra9/F2048²的100 TFF早期诊断值已在主表列出，不能以其Re偏差替代长期统计。',
'','## CPU辅助试验：与RB基准分开解释','',
'这些是已执行的隔离算子/标准解测试，不能称作Ra8或Ra9完整热对流成功。',
'','周期标量波：128节点/周期、均匀周期网格、5000步、指定u=(0.04,0)，分别用对应细/粗级热扩散参数；这里不解浮力耦合，因此Ra仅用于选取参数，不作为本测试物理控制量。误差为测得扩散率对名义扩散率的偏差。',
'','|温度模型|细级扩散偏差%|粗级扩散偏差%|','|---|---:|---:|']
wave=json.loads((C/'Multiblock_algorithm_research_20261020/wave_results.json').read_text())
for a,b in zip(wave[::2],wave[1::2]):lines.append('|{}|{}|{}|'.format(a['model'],fmt(100*a['relative_alpha_error'],5),fmt(100*b['relative_alpha_error'],5)))
lines += ['','同批还做了两级41×41波数扫描、速度幅值0–0.10、两个方向；H-RR温度候选在所扫范围未发现增长模，但没有完成该温度候选的耦合多块RB移植。MRT-A原文率在部分粗级模式不稳定。具体谱半径见算法研究原报告。',
'','流场周期剪切波：均匀周期网格，128节点波长，500–5000步拟合，omega=1.97。名义黏性相对偏差：原算法0.0100%、PR2为0.0195%、RR3为0.0195%、H-RR σ=.98为0.0565%、H-RR σ=.90为0.2108%。这是剪切波误差，不是Nu/Re误差。',
'','双D2Q9导热检查：禁用浮力，N=64，多块F64²或全域均匀64²，推进约5 TFF；Ra=10³仅保留时间标度。修正前温度halo只交换5方向，体积Nu=1、热/冷壁Nu=3.06628/3.16698、中线Nu=.74855；修正为9方向后热/冷壁约.994197、中线约1.000060，与均匀组一致。故早期双D2Q9多块结果有已知缺陷，必须留档但不能作为有效精度依据。',
'','## 复现、探针与程序核验记录','',
'- 原Ra8轨迹从零、检查点续算、CPU/GPU短窗口均复现失稳；逐阶段探针显示上细区粗细交界附近流场短波先放大，后出现非正密度，温度随后损坏。原定期检查约18.508 TFF，带更密探针约18.2136 TFF，频率不同不能混成不同算法收益。',
'- 接口代数400组状态及制造场、148225个空间模板位置、三层时间权重、迁移下标、由f重建rho等检查没有发现所测直接代数错误；不证明全局跨接口守恒与稳定性。',
'- 原诊断Ra8_probe首次编译失败，接口control_a/control_b缺索引的准备失败，UseG首次启动和双D2Q9首次GPU占用退出，均属于运行准备/探针问题，已保留日志；不计成物理发散或成功新算法。',
'- UseG、魔法参数UseG、正则化及双D2Q9分支均有CPU/GPU碰撞和接口对照、小网格推进预检。当前分区χ续算另有N64、Ra8的4 TFF连续/2+2 TFF续算对照：两种χ在CPU、GPU上终态和历史差异为0，错误χ配对被拒绝。这些只能证明所测重启和执行一致性。',
'- 源哈希、参数核对、二进制/最终宏观场有限性记录均保留在各组证据目录。本次新下载的正则化长算逐组重新核对历史与运行源哈希，独立时间积分与原Fortran后处理相符；并读取了终场与完整检查点有限性校验记录。未再次下载所有大型二进制。',
'','## 每一行的说明与出处','']
for r in rows:
    evpath=r['evidence'];full=(REPO/evpath) if evpath.startswith('文档') else (C/evpath)
    lines.append('- **{}** {}{}[证据]({})'.format(r['id'],r['note']+'；' if r['note'] else '', '', full.as_posix()))
lines += ['','辅助报告：',
'- [温度六模型谱/波测试]('+str((C/'Multiblock_algorithm_research_20261020/算法研究与局部测试.md').as_posix())+')',
'- [正则化短测和剪切波]('+str((C/'Multiblock_regularization_20261020/短测结果.md').as_posix())+')',
'- [修正前后导热检查]('+str((C/'Multiblock_LBMCDE_LevelChi_Ra8_Ra9_20261022/conduction_comparison.json').as_posix())+')',
'','## 当前结论与核验边界','',
'**已核查：**本轮算例目录、已有冻结报告、文献原始表格、新读取的服务器长算状态、退出日志、历史和运行源哈希；重新计算正则化长算、加密RR3短测的时间积分与相对误差。',
'','**仍未证明：**多块高Ra的一般长期稳定性、所有方案严格统计收敛、空间网格无关性；仍在运行的Ra9加密RR3最终统计；稳态参考Re定义疑点。',
'','**问题归类：**原始多块与均匀对照差异支持优先研究跨级数值耦合与高频模态；最新分区χ长算还有明显质量守恒缺陷。早期双D2Q9热halo遗漏属于已确认实现问题，已在后续副本修正。不能把全部失败归为后处理，也不能仅由现有数据锁定所有分支的同一根因。',
'','本次仅汇集、读取和复算结果，没有修改求解器，没有新提交或延长算例。']
(P/'全部测试结果.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('Catalog rows:',len(rows))
print('Report:',P/'全部测试结果.md')
