"""Draw the current static block geometry, in fine-grid length units.

PNG and vector SVG are generated from the same geometric primitives.
No solver source is edited. Requires Pillow (available in the bundled runtime).
"""
from pathlib import Path
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / '2DRBOpenaccMultiblock.F90'
src = SOURCE.read_text(encoding='utf-8-sig')
def parameter(name):
    match = re.search(r'\b'+re.escape(name)+r'\s*=\s*(\d+)\b', '\n'.join(line.split('!')[0] for line in src.splitlines()), re.I)
    if not match:
        raise ValueError('Missing integer parameter: '+name)
    return int(match[1])
nx, ny, ratio = [parameter(k) for k in ('nx','ny','refineRatio')]
overlap_cells, fine_ov, skin = [parameter(k) for k in ('coarseOverlapCells','fineOverlapCells','interfaceSkin')]
left,right,bottom,top = [parameter('fineLayerCells'+k) for k in ('Left','Right','Bottom','Top')]
xl,xr,yb,yt = left-.5,nx-right+.5,bottom-.5,ny-top+.5
ov = overlap_cells*ratio
if ratio != 2:
    raise ValueError('The two-substep time diagram requires refineRatio=2; update its stages for another ratio.')
nl,nr,nb,nt=left+fine_ov,right+fine_ov,bottom+fine_ov,top+fine_ov
middle=nx-nl-nr
assert middle>0 and ny-nb-nt>0
assert overlap_cells>=skin and fine_ov>=skin
assert (xr-xl)%ratio==0 and (yt-yb)%ratio==0
blocks=[]
for name,first,counts,dx in [
    ('coarse',[xl-ov,yb-ov],[round((xr-xl)/ratio)+2*overlap_cells+1,
                            round((yt-yb)/ratio)+2*overlap_cells+1],ratio),
    ('left',[.5,.5],[nl,ny],1),
    ('right',[nx-nr+.5,.5],[nr,ny],1),
    ('bottom',[nl+.5,.5],[middle,nb],1),
    ('top',[nl+.5,ny-nt+.5],[middle,nt],1)]:
    last=[first[k]+(counts[k]-1)*dx for k in range(2)]
    blocks.append({'array':name,'dx':dx,'first_computed_node':first,
                   'last_computed_node':last,'computed_nodes':counts,
                   'offset':[x-.5*dx for x in first]})
def fmt(x): return f'{x:g}'

INK='#243147'; MUTED='#536277'; GRID='#d7dee8'
BLUE='#326fc2'; BLUE_FILL='#e8f0fc'; ORANGE='#b96c16'; FINE_FILL='#fff1de'
OVER='#bd3e78'; OVER_FILL='#f5dbe8'; WHITE='#ffffff'
FONT=Path('C:/Windows/Fonts/msyh.ttc')

class Canvas:
    def __init__(self,w,h,title):
        self.im=Image.new('RGB',(w,h),WHITE)
        self.d=ImageDraw.Draw(self.im)
        self.svg=ET.Element('svg',{'xmlns':'http://www.w3.org/2000/svg','width':str(w),'height':str(h),
                                 'viewBox':f'0 0 {w} {h}','role':'img'})
        ET.SubElement(self.svg,'title').text=title
        self.rect(0,0,w,h,WHITE)
    def rect(self,x,y,w,h,fill=None,stroke=None,width=1):
        self.d.rectangle((x,y,x+w,y+h),fill=fill,outline=stroke,width=width)
        ET.SubElement(self.svg,'rect',{'x':str(x),'y':str(y),'width':str(w),'height':str(h),
            'fill':fill or 'none','stroke':stroke or 'none','stroke-width':str(width)})
    def line(self,x1,y1,x2,y2,color=INK,width=2,dash=False):
        if dash:
            dx,dy=x2-x1,y2-y1; length=(dx*dx+dy*dy)**.5
            if length:
                for a in range(0,int(length),15):
                    b=min(a+8,length)
                    self.d.line((x1+dx*a/length,y1+dy*a/length,x1+dx*b/length,y1+dy*b/length),fill=color,width=width)
        else:
            self.d.line((x1,y1,x2,y2),fill=color,width=width)
        attrs={'x1':str(x1),'y1':str(y1),'x2':str(x2),'y2':str(y2),'stroke':color,'stroke-width':str(width)}
        if dash: attrs['stroke-dasharray']='8 7'
        ET.SubElement(self.svg,'line',attrs)
    def text(self,x,y,text,size=30,color=INK,anchor='middle'):
        font=ImageFont.truetype(str(FONT),size)
        for k,line in enumerate(text.split('\n')):
            yy=y+k*(size+10)
            self.d.text((x,yy),line,font=font,fill=color,anchor={'middle':'mm','start':'lm','end':'rm'}[anchor])
            el=ET.SubElement(self.svg,'text',{'x':str(x),'y':str(yy),'font-family':'Microsoft YaHei, Noto Sans CJK SC, sans-serif',
                'font-size':str(size),'fill':color,'text-anchor':anchor,'dominant-baseline':'central'})
            el.text=line
    def dot(self,x,y,r,color,square=False):
        if square:
            self.rect(x-r,y-r,2*r,2*r,color)
        else:
            self.d.ellipse((x-r,y-r,x+r,y+r),fill=color)
            ET.SubElement(self.svg,'circle',{'cx':str(x),'cy':str(y),'r':str(r),'fill':color})
    def arrow(self,x1,y1,x2,y2,color=INK):
        self.line(x1,y1,x2,y2,color,2)
        dx,dy=x2-x1,y2-y1; norm=(dx*dx+dy*dy)**.5; dx/=norm; dy/=norm
        for x,y,s in [(x1,y1,1),(x2,y2,-1)]:
            pts=[(x,y),(x+s*dx*11-dy*5,y+s*dy*11+dx*5),(x+s*dx*11+dy*5,y+s*dy*11-dx*5)]
            self.d.polygon(pts,fill=color)
            ET.SubElement(self.svg,'polygon',{'points':' '.join(f'{a},{b}' for a,b in pts),'fill':color})
    def save(self,name):
        self.im.save(HERE/(name+'.png'))
        ET.indent(self.svg)
        ET.ElementTree(self.svg).write(HERE/(name+'.svg'),encoding='utf-8',xml_declaration=True)

SAME='#237b69'
layout=Canvas(1500,1500,'中心粗网格与四套紧凑细数组；粗细两侧按各自格距延伸。')
layout.text(750,48,'五套数组：粗细网格采用非对称重叠',44)
layout.text(750,105,f'nx={nx}，ny={ny}  ·  粗细比 {ratio}  ·  坐标以细格距计',29,MUTED)
for x,color,desc in [(175,BLUE_FILL,'中心粗网格'),(620,FINE_FILL,'四套细数组'),(1040,OVER_FILL,'粗细重叠带')]:
    layout.rect(x,150,28,24,color,GRID); layout.text(x+43,162,desc,27,anchor='start')
X0,Y0=300,245
sx=sy=min(900/nx,900/ny)
def X(x): return X0+x*sx
def Y(y): return Y0+(ny-y)*sy
def box(xlo,xhi,ylo,yhi,fill): layout.rect(X(xlo),Y(yhi),(xhi-xlo)*sx,(yhi-ylo)*sy,fill)
box(0,nx,0,ny,FINE_FILL);box(xl,xr,yb,yt,BLUE_FILL)
# Intersection of the coarse rectangle and fine ring, using node coordinates.
box(xl-ov,xr+ov,yb-ov,yb+fine_ov,OVER_FILL)
box(xl-ov,xr+ov,yt-fine_ov,yt+ov,OVER_FILL)
box(xl-ov,xl+fine_ov,yb+fine_ov,yt-fine_ov,OVER_FILL)
box(xr-fine_ov,xr+ov,yb+fine_ov,yt-fine_ov,OVER_FILL)
layout.rect(X(0),Y(ny),nx*sx,ny*sy,None,INK,3)
for y in (yb,yt):
    layout.line(X(xl),Y(y),X(xr),Y(y),OVER,3,True)
for x in (xl,xr): layout.line(X(x),Y(yb),X(x),Y(yt),OVER,3,True)
# Fine-array storage seams are halfway between adjacent fine nodes.
for x in (nl,nx-nr):
    for ya,ye in [(0,nb),(ny-nt,ny)]: layout.line(X(x),Y(ya),X(x),Y(ye),SAME,3)
layout.text(X(nx/2),Y((yt+ny)/2),f'top：{middle} × {nt}，dx=1',29)
layout.text(X(nx/2),Y(yb/2),f'bottom：{middle} × {nb}，dx=1',29)
layout.text(X(xl/2),Y(ny*.57),f'left\n{nl} ×\n{ny}\ndx=1',25)
layout.text(X((xr+nx)/2),Y(ny*.57),f'right\n{nr} ×\n{ny}\ndx=1',25)
layout.text(X((xl+xr)/2),Y(ny*.64),'coarse：中心粗网格',40)
layout.text(X((xl+xr)/2),Y(ny*.56),f'dx={ratio}，Δt={ratio}',32)
layout.text(X((xl+xr)/2),Y(ny*.49),f'积分范围 [{fmt(xl)}, {fmt(xr)}] × [{fmt(yb)}, {fmt(yt)}]',26)
layout.text(X((xl+xr)/2),Y(ny*.42),f'物理宽高 {fmt(xr-xl)} × {fmt(yt-yb)}',29)
layout.text(X((xl+xr)/2),Y(ny*.34),f'计算节点 {blocks[0]["computed_nodes"][0]} × {blocks[0]["computed_nodes"][1]}（含重叠）',28)
layout.text(X((xl+xr)/2),Y(ny*.27),f'粗向外延伸 {overlap_cells} 个粗格距；细向内延伸 {fine_ov} 个细格距',25,MUTED)
for tick in (0,xl,xr,nx): layout.text(X(tick),1160,fmt(tick),25)
for tick in (0,yb,yt,ny): layout.text(275,Y(tick),fmt(tick),25,anchor='end')
layout.text(1305,1160,'x',27);layout.text(275,212,'y',27)
layout.line(100,1220,155,1220,OVER,4,True);layout.text(178,1220,'粉色带：粗细共同计算范围；虚线：积分分区线',29,anchor='start')
layout.line(100,1275,155,1275,SAME,4)
layout.text(178,1275,'绿色短线：细数组接缝，直接交换 f_post / g_post 后迁移',29,SAME,anchor='start')
layout.text(750,1330,f'Left={left} → x={fmt(xl)}     Right={right} → x={fmt(xr)}',28)
layout.text(750,1378,f'Bottom={bottom} → y={fmt(yb)}     Top={top} → y={fmt(yt)}',28)
layout.text(750,1443,'Left / Right / Bottom / Top 是从对应墙面数的节点编号，不是物理层厚。',26,MUTED)
layout.save('multiblock-layout')

detail=Canvas(1500,1360,'左侧非对称重叠：区分积分分区线、人工边缘和两层接收节点。')
detail.text(750,45,f'左侧接口放大：粗延伸 {overlap_cells} 格，细延伸 {fine_ov} 格',42)
lo,hi=xl-ov,xl+fine_ov
xa,xe=xl-14,xl+18
ylo=blocks[0]['first_computed_node'][1]+ratio*round((ny/2-6-blocks[0]['first_computed_node'][1])/ratio)
yhi=ylo+12
scale=34
xx=lambda x: 190+(x-xa)*scale
yy=lambda y: 280+(yhi-y)*scale
detail.text(750,100,f'左细块与中心粗块；分界 x={fmt(xl)}，粗细两块在同坐标各存一份状态',27,MUTED)
detail.rect(xx(xa),yy(yhi),(xl-xa)*scale,12*scale,FINE_FILL)
detail.rect(xx(xl),yy(yhi),(xe-xl)*scale,12*scale,BLUE_FILL)
detail.rect(xx(lo),yy(yhi),(hi-lo)*scale,12*scale,OVER_FILL)
for i in range(round(hi-xa)+1):
    x=xa+i
    detail.line(xx(x),yy(ylo),xx(x),yy(yhi),'#e4cbb8',1)
    for j in range(13): detail.dot(xx(x),yy(ylo+j),4,ORANGE)
for i in range(round((xe-lo)/ratio)+1):
    x=lo+i*ratio
    detail.line(xx(x),yy(ylo),xx(x),yy(yhi),'#a9bfdf',1)
    for j in range(0,13,ratio): detail.rect(xx(x)-8,yy(ylo+j)-8,16,16,None,BLUE,2)
for j in range(13): detail.line(xx(xa),yy(ylo+j),xx(hi),yy(ylo+j),'#e4cbb8',1)
for j in range(0,13,ratio): detail.line(xx(lo),yy(ylo+j),xx(xe),yy(ylo+j),'#a9bfdf',1)
detail.line(xx(xl),yy(ylo)-5,xx(xl),yy(yhi)-5,OVER,3,True)
detail.arrow(xx(lo),190,xx(hi),190,OVER)
detail.text(750,150,f'重叠范围 [{fmt(lo)}, {fmt(hi)}]，跨度 {ov+fine_ov} 个细格距',28,OVER)
detail.text(xx(lo)-18,242,f'粗人工边缘 {fmt(lo)}',25,anchor='end')
detail.text(xx(hi)+18,242,f'细人工边缘 {fmt(hi)}',25,anchor='start')
for x in (xa,lo,xl,hi,xe):
    tick_y=758 if x==hi else 728
    detail.line(xx(x),yy(ylo)+12,xx(x),tick_y-18,MUTED,1)
    detail.text(xx(x),tick_y,fmt(x),25)
for y in (ylo,ylo+6,yhi): detail.text(165,yy(y),fmt(y),25,anchor='end')
detail.dot(200,795,4,ORANGE);detail.text(225,795,'细点 dx=1',27,anchor='start')
detail.rect(560,787,16,16,None,BLUE,2);detail.text(590,795,f'粗点 dx={ratio}',27,anchor='start')
detail.rect(980,787,16,16,None,BLUE,2);detail.dot(988,795,4,ORANGE);detail.text(1010,795,'粗细重合点',27,anchor='start')
detail.text(750,860,f'虚线 x={fmt(xl)}：积分分区线；两侧人工边缘仍在流体内部，不是墙壁',27)
detail.text(750,915,f'粗向左：{overlap_cells} × {ratio} = {ov} 细格距；细向右：{fine_ov} 细格距',28)
detail.rect(110,959,1280,260,'#f5f7fa',GRID)
detail.text(140,997,f'interfaceSkin={skin}：各自接收 {skin} 层，重建 f、g 后参加碰撞与迁移',28,anchor='start')
cs='、'.join(fmt(lo+i*ratio) for i in range(skin));fs='、'.join(fmt(hi-skin+1+i) for i in range(skin))
detail.text(140,1055,f'粗块接收列：{cs}     |     细块接收列：{fs}',27,anchor='start')
for i in range(skin):
    for j in range(0,13,ratio): detail.rect(xx(lo+i*ratio)-9,yy(ylo+j)-9,18,18,None,OVER,3)
    for j in range(13): detail.dot(xx(hi-i),yy(ylo+j),6,SAME)
detail.text(140,1110,'图中粉框：粗接收节点；绿点：细接收节点（其余为正常计算节点）',26,anchor='start')
detail.text(140,1165,'共址取值也要转换非平衡矩尺度；接收层不作为反向交换的来源。',26,anchor='start')
detail.text(750,1265,'空间：共址直接取值，非共址每方向四点 Lagrange；模板可伸出重叠区。',26,MUTED)
detail.text(750,1310,'所有坐标以细格距计；人工边缘表示最外计算节点，没有额外半格距壁面。',26,MUTED)
detail.save('multiblock-interface')

same=Canvas(1500,820,'四套细数组通过碰撞后分布函数外圈交换，保持同级连续迁移。')
same.text(750,48,'四套细数组：接缝处先交换，再迁移',42)
same.text(750,107,'left / right / bottom / top 各存一套 f、g；格距和时间步相同',28,MUTED)
same.rect(100,185,1300,190,FINE_FILL,GRID)
same.text(750,240,'先将邻区 f_post / g_post 复制到本区迁移外圈',32)
same.text(750,310,'f(i,j,α) ← f_post(i−ex(α), j−ey(α), α)',31)
same.text(750,440,'同级接缝不插值、不反弹；D2Q9 同时交换对角方向所需数据。',29,SAME)
same.text(750,515,'顺序：四区流场碰撞 → 交换 → 迁移；温度场随后按同样顺序推进。',28)
same.text(750,605,'左右细区贯穿全高，上下细区填中间；角点不重复存储。',29)
same.text(750,715,'中心空区不分配细数组；粗细人工边缘仍用 interfaceSkin 接收层交换。',27,MUTED)
same.save('multiblock-samelevel')

clock=Canvas(1500,1170,'一次粗步：先准备碰撞前缓冲，粗步为2、细步为1，在t+2同步。')
clock.text(750,48,'读懂一次粗步：数据先到，节点再碰撞',43)
clock.text(750,104,'refineRatio = 2    |    以下时间以一个细步为单位',29,MUTED)
stages=[
    (170,'入口','两级都在 t；初始交换或上次同步已准备好缓冲',
     '交换量来自迁移后的 f、g，供接下来的碰撞使用', '#f1f4f8'),
    (345,'① 粗块预测','粗块：碰撞 → 迁移 → 原温度推进，到达 t+2',
     '保留 t−2、t、t+2 三个粗时间层；来源点排除两层人工边界',BLUE_FILL),
    (520,'② 细块第一步','准备 t 的接口 → 碰撞、迁移及温度推进 → t+1',
     '重建细接收层；四套细数组先交换碰撞后外圈，再迁移',FINE_FILL),
    (695,'③ 细块第二步','准备 t+1 的接口 → 碰撞、迁移及温度推进 → t+2',
     '粗历史插值到 t+1；细数组间继续交换碰撞后外圈',FINE_FILL),
    (870,'④ 同步','细 → 粗修复粗缓冲；再补齐细缓冲，均在 t+2',
     '滚动历史，时钟增加2；此时输出、存检查点或开始下一个粗步',OVER_FILL),
]
for y,title,action,note,fill in stages:
    clock.rect(55,y,1390,137,fill,GRID)
    clock.text(82,y+37,title,29,anchor='start')
    clock.text(330,y+38,action,29,anchor='start')
    clock.text(330,y+91,note,25,MUTED,anchor='start')
    if y<870: clock.text(750,y+153,'↓',29,MUTED)
clock.text(750,1070,'中间时刻：Q(t+1) = −Q(t−2)/8 + 3Q(t)/4 + 3Q(t+2)/8',29)
clock.text(750,1125,'首个粗步缺少 t−2，使用 t 与 t+2 的平均值启动；收到的缓冲节点必须参加碰撞。',26,MUTED)
clock.save('multiblock-timestep')

manifest={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'nx':nx,'ny':ny,
          'refineRatio':ratio,'fineLayerCellsLeft':left,'fineLayerCellsRight':right,
          'fineLayerCellsBottom':bottom,'fineLayerCellsTop':top,'coarseOverlapCells':overlap_cells,
          'fineOverlapCells':fine_ov,
          'interfaces':{'xLeft':xl,'xRight':xr,'yBottom':yb,'yTop':yt},
          'coarse_extension_in_fine_units':ov,'fine_extension_in_fine_units':fine_ov,
          'shared_node_span':ov+fine_ov,'interfaceSkin':skin,'blocks':blocks,
          'computed_node_total':sum(b['computed_nodes'][0]*b['computed_nodes'][1] for b in blocks),
          'integration_sample_total':None,
          'legend':{'pink':'coarse-fine overlap; dashed line is integration partition',
                    'fine_seams':'green; exchange post-collision halos then stream'}}
for b in blocks:
    area=0.; positive=0; dx=b['dx']
    for j in range(b['computed_nodes'][1]):
        y=b['first_computed_node'][1]+j*dx
        for i in range(b['computed_nodes'][0]):
            x=b['first_computed_node'][0]+i*dx
            cut=max(0,min(x+dx/2,xr)-max(x-dx/2,xl))*max(0,min(y+dx/2,yt)-max(y-dx/2,yb))
            weight=cut if b['array']=='coarse' else 1-cut
            assert weight>=0
            area+=weight; positive+=int(weight>0)
    b['integration_area']=area
    b['positive_area_node_count']=positive
assert sum(b['integration_area'] for b in blocks)==nx*ny
manifest['active_node_total']=manifest['computed_node_total']
manifest['integration_sample_total']=sum(b['positive_area_node_count'] for b in blocks)
manifest['fine_node_total']=sum(b['computed_nodes'][0]*b['computed_nodes'][1] for b in blocks[1:])
(HERE/'layout_parameters.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
