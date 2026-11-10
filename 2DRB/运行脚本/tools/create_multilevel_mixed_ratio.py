"""Derive the mixed-interface-ratio solver without modifying its approved parent."""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[2]
PARENT=ROOT/'多块网格/2DRBOpenaccMultilevel.F90'
TARGET=ROOT/'多块网格/2DRBOpenaccMultilevelMixedRatio.F90'
OUT=ROOT/'运行脚本/run_cases/Multilevel_mixed_ratio_implementation_20261030'
raw=PARENT.read_bytes()
parent_sha=hashlib.sha256(raw).hexdigest()
assert parent_sha=='0c24f5e544ff1308587f7ec7381d4a7e9066b995255424a5ab525f6636672ce5'
text=raw.decode('utf-8-sig').replace('\r\n','\n')
changes=[]
def replace(old,new,count=1):
    global text
    assert text.count(old)==count,(old,text.count(old),count)
    text=text.replace(old,new)
    changes.append(dict(old=old,new=new))

replace('!!! 来源：2DRBOpenaccMultiblock.F90 的独立副本；原文件保持不变。',
        '!!! 来源：2DRBOpenaccMultilevel.F90 的独立混合比例副本；原文件保持不变。')
replace('!!! 原文件 SHA256：84ca3661eaa2544cca9bb19403b1eaa7d97e82aa4f09ccbbd7cf657de915aa2c',
        '!!! 父本 SHA256：'+parent_sha)
replace('!!! 第1级贴壁且最细；向中心每增加一级，dx、dt均乘2。',
        '!!! 第1级贴壁且最细；向中心按各接口比例增大dx、dt，默认4倍再2倍。')
replace('        ! 默认三级：壁面1024级 -> 中间512级 -> 中心256级。\n'
        '        ! 只改nGridLevels即可切换为2级或4级；1级是均匀网格回归分支。\n'
        '        ! nx、ny以最细格距计，必须能被2**(nGridLevels-1)整除。',
        '        ! 默认三级：壁面4096级 -> 中间1024级 -> 中心512级。\n'
        '        ! nGridLevels包含中心级；改变级数时，同时填写每个接口的比例数组。\n'
        '        ! 例：两级[4]，三级[4,2]或[2,2]，四级[2,2,2]；一级填[1]作为占位。\n'
        '        ! nx、ny以最细格距计，必须能被全部接口比例的乘积整除。')
replace('integer(kind=4), parameter :: nx=1024, ny=1024','integer(kind=4), parameter :: nx=4096, ny=4096')
replace('        integer(kind=4), parameter :: levelRefineRatio=2    ! 相邻级固定为2倍，避免直接跨越多级\n'
        '        ! 固定相邻比2，位移等价于2**max(0,nGridLevels-1)，绕过NVHPC 24.3的幂常量折叠问题。\n'
        '        integer(kind=4), parameter :: coarsestStep=ishft(1,max(0,nGridLevels-1))',
        '        integer(kind=4), parameter :: interfaceRefineRatio(max(1,nGridLevels-1))=[4,2]\n'
        '        ! 比例按最细向中心排列，仅接受2或4；累计格距/时间步默认为1、4、8。\n'
        '        ! 一级占位1避免NVHPC 24.3归约零长常量数组；外层max让非法非正比例到初始化时报错。\n'
        '        integer(kind=4), parameter :: coarsestStep=max(1,product(interfaceRefineRatio))')
replace('! 第1级的环带厚度按最细格子数计；64层细格作为初始选择，不代表已完成网格独立验证。',
        '! 第1级环带按最细格子计；默认256格，厚度为边长1/16，尚未完成网格独立验证。')
replace('finestLayerCellsLeftRight=64','finestLayerCellsLeftRight=256')
replace('finestLayerCellsBottomTop=64','finestLayerCellsBottomTop=256')
replace('! 沿用原共址界面约定：默认左右分区线为63.5与959.5，两侧厚度相差一个最细格距。',
        '! 沿用原共址界面约定：默认左右分区线为255.5与3839.5，两侧相差一个最细格距。')
replace('! 默认第2级32格对应64个最细格距；第3级32格对应128个最细格距。',
        '! 默认第2级64格乘本级格距4，对应256个最细格距；中心级不使用环带厚度。')
replace('innerLayerCellsLeftRight=32','innerLayerCellsLeftRight=64')
replace('innerLayerCellsBottomTop=32','innerLayerCellsBottomTop=64')
replace('Rayleigh=1.0d7','Rayleigh=1.0d9')
replace("settingsFile='SimulationSettings2DOpenaccMultilevel.txt'","settingsFile='SimulationSettings2DOpenaccMultilevelMixedRatio.txt'")
replace("reloadFilePrefix='reloadFile2DOpenaccMultilevel'","reloadFilePrefix='reloadFile2DOpenaccMultilevelMixedRatio'")
replace("        if(levelRefineRatio/=2) error stop 'This nested hierarchy uses adjacent refinement ratio two'",
        "        if(nGridLevels==1) then\n"
        "            if(any(interfaceRefineRatio/=1)) error stop 'A single level requires the ratio placeholder one'\n"
        "        else\n"
        "            if(any(interfaceRefineRatio/=2.and.interfaceRefineRatio/=4)) &\n"
        "                error stop 'Each interface refinement ratio must be two or four'\n"
        "        endif")
replace("error stop 'nx,ny must be divisible by 2**(nGridLevels-1)'",
        "error stop 'nx,ny must be divisible by the product of interface refinement ratios'")
replace('            dx=dble(levelRefineRatio**(l-1))',
        '            dx=1.0d0\n            if(l>1) dx=grid(l-1)%dx*dble(interfaceRefineRatio(l-1))')
replace('                if(min(tx,ty)<2*interfaceSkin+4) &\n'
        "                    error stop 'Ring thickness must be at least 2*interfaceSkin+4 local cells'",
        '                if(min(tx,ty)<interfaceSkin*interfaceRefineRatio(l)+4) &\n'
        "                    error stop 'Ring thickness must be at least interfaceSkin*localRatio+4 local cells'")
replace('if(min(ix2-ix1,iy2-iy1)<levelRefineRatio*3)',
        'if(min(ix2-ix1,iy2-iy1)<interfaceRefineRatio(l)*3)')
replace('if(mod(ix2-ix1,levelRefineRatio)/=0.or.mod(iy2-iy1,levelRefineRatio)/=0)',
        'if(mod(ix2-ix1,interfaceRefineRatio(l))/=0.or.mod(iy2-iy1,interfaceRefineRatio(l))/=0)')
replace('! 递归同步：先预测本级终点，再让外侧更细级走两个子步，最后双向补齐接收层。',
        '! 递归同步：先预测本级终点，再按本接口比例走2或4个细子步，最后双向补齐接收层。')
replace('        do k=1,levelRefineRatio\n            theta=dble(k-1)/dble(levelRefineRatio)',
        '        do k=1,interfaceRefineRatio(l-1)\n'
        '            theta=dble(k-1)/dble(interfaceRefineRatio(l-1))')
replace("ML2DSNAPSHOT0001","ML2DSNAPSHOT0002")
replace('write(k) magic,current_config(),current_physics(),reloadItcOffset+itc,nBlocks,timeUnit',
        'write(k) magic,current_config(),current_physics(),interfaceRefineRatio,reloadItcOffset+itc,nBlocks,timeUnit')
replace("write(k,'(a)') 'Parent SHA256: 84ca3661eaa2544cca9bb19403b1eaa7d97e82aa4f09ccbbd7cf657de915aa2c'",
        "write(k,'(a)') 'Parent SHA256: "+parent_sha+"'")
replace("        write(k,'(a,16(i0,1x))') 'Configuration: ',current_config()",
        "        write(k,'(a,16(i0,1x))') 'Configuration: ',current_config()\n"
        "        write(k,'(a,*(i0,1x))') 'Interface refinement ratios (fine to coarse): ',interfaceRefineRatio")
replace('ML2DRESTART0001','ML2DRESTART0002',count=4)
replace('write(k) magic,current_config(),current_physics(),total,NuReSampleCount,errorU,errorT,nBlocks',
        'write(k) magic,current_config(),current_physics(),interfaceRefineRatio,total,NuReSampleCount,errorU,errorT,nBlocks')
replace('integer :: k,b,l,ios,config(16),total,samples,nb,geometry(5),nf,cycles',
        'integer :: k,b,l,ios,config(16),ratios(max(1,nGridLevels-1)),total,samples,nb,geometry(5),nf,cycles')
replace('read(k,iostat=ios) magic,config,physics,total,samples,savedU,savedT,nb',
        'read(k,iostat=ios) magic,config,physics\n'
        "        if(ios/=0) error stop 'Truncated multilevel checkpoint configuration'\n"
        "        if(magic/='ML2DRESTART0002') error stop 'Checkpoint is not the mixed-ratio multilevel format'\n"
        "        if(any(config/=current_config())) error stop 'Checkpoint mesh/BC configuration mismatch'\n"
        '        read(k,iostat=ios) ratios,total,samples,savedU,savedT,nb')
replace("        if(any(config/=current_config()).or.nb/=nBlocks) error stop 'Checkpoint mesh/BC configuration mismatch'",
        "        if(any(ratios/=interfaceRefineRatio).or.nb/=nBlocks) error stop 'Checkpoint interface ratio configuration mismatch'")
replace('cycles/=total/levelRefineRatio**l','cycles/=total/nint(grid(l+1)%dx)')
replace('! 多级时dt与四个松弛率均取当前来源级；第1级dt=1，第2级dt=2，依此类推。',
        '! dt与四个松弛率均取当前来源级；默认[4,2]时dt=1、4、8，全2比例时为1、2、4。')
replace('! 以第2级来源dt=2为例，括号内为m_a-meq_a+(F_m,a)；细块dt=1时为m_a-meq_a+(F_m,a)/2。',
        '! 全2比例的第2级dt=2示例中，括号内为m_a-meq_a+(F_m,a)；默认第2级dt=4时为m_a-meq_a+2*(F_m,a)。')
replace('! 第2级来源的缩放系数为s(a)/2；buffer(6)、buffer(9)、buffer(11)对应守恒矩，有限场下为零。',
        '! 来源级缩放系数为s(a)/dt；buffer(6)、buffer(9)、buffer(11)对应守恒矩，有限场下为零。')
replace('! 8. buffer(15:19)保存5个缩放后的温度非平衡矩；第2级来源缩放系数为q(a)/2；一般级为q(a)/dt。',
        '! 8. buffer(15:19)保存5个缩放后的温度非平衡矩；来源级缩放系数为q(a)/dt。')
replace('! 4. 用接收网格的松弛率和时间步恢复完整矩；第1级到第2级时dt=2，第2级到第3级时dt=4。',
        '! 4. 用接收级的松弛率和时间步恢复完整矩；默认[4,2]的第2级dt=4，第3级dt=8。')
assert 'levelRefineRatio' not in text
assert not TARGET.exists(),'Do not overwrite a developed or submitted solver'
TARGET.write_text(text,encoding='utf-8',newline='\n')
OUT.mkdir(exist_ok=True)
(OUT/'parent.F90').write_bytes(raw)
(OUT/'derivation.json').write_text(json.dumps(dict(parent=str(PARENT),target=str(TARGET),parent_sha256=parent_sha,target_sha256=hashlib.sha256(TARGET.read_bytes()).hexdigest(),changes=changes),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(source=str(TARGET),changes=len(changes),sha256=hashlib.sha256(TARGET.read_bytes()).hexdigest())))
