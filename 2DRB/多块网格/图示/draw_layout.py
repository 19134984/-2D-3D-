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
overlap_cells, skin = [parameter(k) for k in ('overlapCells','interfaceSkin')]
left,right,bottom,top = [parameter('fineLayerCells'+k) for k in ('Left','Right','Bottom','Top')]
xl,xr,yb,yt = left-.5,nx-right+.5,bottom-.5,ny-top+.5
ov = overlap_cells*ratio
if ratio != 2:
    raise ValueError('The two-substep time diagram requires refineRatio=2; update its stages for another ratio.')
blocks=[]
for name,xa,xb,ya,ye,h in [('1 中心',xl,xr,yb,yt,ratio),('2 连通细环',0,nx,0,ny,1)]:
    first=[max(.5,xa-ov),max(.5,ya-ov)]
    last=[min(nx-.5,xb+ov),min(ny-.5,ye+ov)]
    counts=[round((last[k]-first[k])/h)+1 for k in range(2)]
    owned_counts=[]
    for axis,(lo,hi) in enumerate([(xa,xb),(ya,ye)]):
        assert abs(first[axis]+(counts[axis]-1)*h-last[axis])<1e-10
        weights=[max(0,min(hi,first[axis]+i*h+h/2)-max(lo,first[axis]+i*h-h/2)) for i in range(counts[axis])]
        assert abs(sum(weights)-(hi-lo))<1e-10
        owned_counts.append(sum(w>0 for w in weights))
    blocks.append({'block':name,'owned_faces':[xa,xb,ya,ye],'h':h,
                   'first_computed_node':first,'last_computed_node':last,
                   'integration_nodes':owned_counts,'computed_nodes':counts})
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
layout=Canvas(1500,1500,'中心粗网格和连通细网格环：仅在两级之间保留接口。')
layout.text(750,48,'连通细网格环：只保留粗细接口',44)
layout.text(750,105,f'nx={nx}，ny={ny}  ·  粗细比 {ratio}  ·  坐标以细格距计',29,MUTED)
for x,color,desc in [(175,BLUE_FILL,'中心粗网格'),(620,FINE_FILL,'一个连通细环'),(1040,OVER_FILL,'粗细重叠带')]:
    layout.rect(x,150,28,24,color,GRID); layout.text(x+43,162,desc,27,anchor='start')
X0,Y0=300,245
sx=sy=min(900/nx,900/ny)
def X(x): return X0+x*sx
def Y(y): return Y0+(ny-y)*sy
def box(xlo,xhi,ylo,yhi,fill): layout.rect(X(xlo),Y(yhi),(xhi-xlo)*sx,(yhi-ylo)*sy,fill)
box(0,nx,0,ny,FINE_FILL);box(xl,xr,yb,yt,BLUE_FILL)
# Coarse/fine overlap occupies only the four sides of the central block.
# Fine/fine buffer overlaps are identified by separate green seams, not this pink category.
for y in (yb,yt): box(xl-ov,xr+ov,y-ov,y+ov,OVER_FILL)
for x in (xl,xr): box(x-ov,x+ov,yb-ov,yt+ov,OVER_FILL)
layout.rect(X(0),Y(ny),nx*sx,ny*sy,None,INK,3)
for y in (yb,yt):
    layout.line(X(xl),Y(y),X(xr),Y(y),OVER,3)
    # 原细块接缝已消除；这里没有人工边界或交换线。
for x in (xl,xr): layout.line(X(x),Y(yb),X(x),Y(yt),OVER,3)
layout.text(X(nx/2),Y((yt+ny)/2),'2  连通细网格环   h=1，Δt=1',30)
layout.text(X(nx/2),Y(yb/2),'同一套 f、g，沿细环直接迁移',30)
layout.text(X(xl/2),Y(ny*.57),'细环\nh=1',27)
layout.text(X((xr+nx)/2),Y(ny*.57),'细环\nh=1',27)
layout.text(X((xl+xr)/2),Y(ny*.64),'1  中心粗块',43)
layout.text(X((xl+xr)/2),Y(ny*.56),f'h={ratio}，Δt={ratio}',32)
layout.text(X((xl+xr)/2),Y(ny*.49),f'基准范围 [{fmt(xl)}, {fmt(xr)}] × [{fmt(yb)}, {fmt(yt)}]',26)
layout.text(X((xl+xr)/2),Y(ny*.42),f'物理宽高 {fmt(xr-xl)} × {fmt(yt-yb)}',29)
layout.text(X((xl+xr)/2),Y(ny*.34),f'计算节点 {blocks[0]["computed_nodes"][0]} × {blocks[0]["computed_nodes"][1]}（含重叠）',28)
layout.text(X((xl+xr)/2),Y(ny*.27),f'正面积积分节点 {blocks[0]["integration_nodes"][0]} × {blocks[0]["integration_nodes"][1]}',27,MUTED)
for tick in (0,xl,xr,nx): layout.text(X(tick),1160,fmt(tick),25)
for tick in (0,yb,yt,ny): layout.text(275,Y(tick),fmt(tick),25,anchor='end')
layout.text(1305,1160,'x',27);layout.text(275,212,'y',27)
layout.line(100,1220,155,1220,OVER,4);layout.text(178,1220,'粉色：粗细接口及重叠带，需要跨级取值与矩重标定',29,anchor='start')
layout.text(178,1275,'四角处连续相邻：没有细块拼接线，没有同级缓冲交换',29,SAME,anchor='start')
layout.text(750,1330,f'Left={left} → x={fmt(xl)}     Right={right} → x={fmt(xr)}',28)
layout.text(750,1378,f'Bottom={bottom} → y={fmt(yb)}     Top={top} → y={fmt(yt)}',28)
layout.text(750,1443,'Left / Right / Bottom / Top 是从对应墙面数的节点编号，不是物理层厚。',26,MUTED)
layout.save('multiblock-layout')

detail=Canvas(1500,1200,'左细块与中心粗块：基准边界就是共址节点列，图中展示实际格点。')
detail.text(750,45,'粗细接口放大：基准边界就是共同节点列',42)
lo,hi=xl-ov,xl+ov
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
detail.text(750,150,f'重叠范围 [{fmt(lo)}, {fmt(hi)}]，跨度 {2*ov}：{2*ov+1} 排细点、{2*overlap_cells+1} 排粗点',28,OVER)
detail.text(xx(lo)-18,242,f'粗块首列 {fmt(lo)}',25,anchor='end')
detail.text(xx(hi)+18,242,f'细块末列 {fmt(hi)}',25,anchor='start')
for x in (xa,lo,xl,hi,xe): detail.text(xx(x),728,fmt(x),25)
for y in (ylo,ylo+6,yhi): detail.text(165,yy(y),fmt(y),25,anchor='end')
detail.dot(200,795,4,ORANGE);detail.text(225,795,'细点 h=1',27,anchor='start')
detail.rect(560,787,16,16,None,BLUE,2);detail.text(590,795,f'粗点 h={ratio}',27,anchor='start')
detail.rect(980,787,16,16,None,BLUE,2);detail.dot(988,795,4,ORANGE);detail.text(1010,795,'粗细重合点',27,anchor='start')
detail.text(750,860,f'虚线 x={fmt(xl)}：共用基准边界，同时也是积分分界',28)
detail.text(750,915,f'每块向外延伸 overlapCells × refineRatio = {overlap_cells} × {ratio} = {ov}',28)
detail.rect(110,959,1280,146,'#f5f7fa',GRID)
detail.text(140,997,f'interfaceSkin={skin}：人工边界接收 {skin} 层，接收后参加碰撞',28,anchor='start')
cs='、'.join(fmt(lo+i*ratio) for i in range(skin));fs='、'.join(fmt(hi-skin+1+i) for i in range(skin))
detail.text(140,1055,f'粗块接收列：{cs}     |     细块接收列：{fs}',27,anchor='start')
detail.text(750,1150,'重合点直接读取交换量；非重合点四点 Lagrange（二维 4×4），随后重建分布。',26,MUTED)
detail.save('multiblock-interface')

same=Canvas(1500,820,'原上细块和左细块现已属于同一细环，直接访问相邻节点。')
same.text(750,48,'原来的上、左细区：现在直接连续迁移',42)
same.text(750,107,'只有一套细网格 f、g；没有同级连接、缓冲打包或插值',28,MUTED)
same.rect(100,185,1300,190,FINE_FILL,GRID)
same.text(750,240,'同一个细网格数组内的普通相邻节点',32)
same.text(750,310,'f(i,j,α) ← f_post(i−ex(α), j−ey(α), α)',31)
same.text(750,440,'穿过原拼接位置与其他细网格内部位置完全相同。',31,SAME)
same.text(750,515,'仍保留粗细交换：中心粗网格 ↔ 连通细环的内边缘。'.replace('↔','与'),29)
same.text(750,605,'二维细数组覆盖全域；中心空区不碰撞、不迁移、不计入积分。',29)
same.text(750,715,'数组为空区保留存储；这版尚未采用只存环上节点的紧凑布局。',27,MUTED)
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
     '细环内边缘读取粗块 t 快照；环内邻点直接迁移',FINE_FILL),
    (695,'③ 细块第二步','准备 t+1 的接口 → 碰撞、迁移及温度推进 → t+2',
     '细环内边缘由粗时间层插值到 t+1；没有同级交换',FINE_FILL),
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
          'fineLayerCellsBottom':bottom,'fineLayerCellsTop':top,'overlapCells':overlap_cells,
          'interfaces':{'xLeft':xl,'xRight':xr,'yBottom':yb,'yTop':yt},
          'extension_from_reference_node':ov,'shared_node_span':2*ov,'interfaceSkin':skin,'blocks':blocks,
          'computed_node_total':sum(b['computed_nodes'][0]*b['computed_nodes'][1] for b in blocks),
          'integration_sample_total':None,
          'legend':{'pink':'coarse-fine interface and overlap','fine_seams':'none; direct neighbor streaming within one ring'}}
fine_active_count=0;fine_owned_count=0
for j in range(ny):
    y=j+.5
    for i in range(nx):
        x=i+.5
        fine_active_count+=int(x<=xl+ov or x>=xr-ov or y<=yb+ov or y>=yt-ov)
        cut=max(0,min(x+.5,xr)-max(x-.5,xl))*max(0,min(y+.5,yt)-max(y-.5,yb))
        fine_owned_count+=int(1-cut>0)
manifest['active_node_total']=blocks[0]['computed_nodes'][0]*blocks[0]['computed_nodes'][1]+fine_active_count
manifest['integration_sample_total']=blocks[0]['integration_nodes'][0]*blocks[0]['integration_nodes'][1]+fine_owned_count
manifest['blocks'][1]['active_node_count']=fine_active_count
manifest['blocks'][1]['positive_area_node_count']=fine_owned_count
manifest['blocks'][1]['owned_faces_note']='Bounding rectangle only; subtract the central owned rectangle in 2D.'
(HERE/'layout_parameters.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
