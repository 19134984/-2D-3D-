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
nx, ny = map(int, re.search(r'parameter :: nx=(\d+), ny=(\d+)', src).groups())
ratio = int(re.search(r'parameter :: refineRatio=(\d+)', src)[1])
overlap_cells = int(re.search(r'parameter :: overlapCells=(\d+)', src)[1])
skin = int(re.search(r'parameter :: interfaceSkin=(\d+)', src)[1])
wx, wy = map(int, re.search(r'parameter :: fineLayerCellsLeftRight=(\d+), fineLayerCellsBottomTop=(\d+)', src).groups())
ov = overlap_cells * ratio
assert (nx, ny, ratio, wx, wy, ov, skin) == (1024, 1024, 2, 128, 128, 4, 2), \
    'This annotated composition targets the current defaults; update labels if the source geometry changes.'

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

layout=Canvas(1500,1440,'当前五块网格：四周细网格、中心粗网格；计算区跨分区接口相互重叠。')
layout.text(750,48,'当前多块网格布局',46)
layout.text(750,105,'nx = ny = 1024    ·    refineRatio = 2    ·    坐标以最细格距计',29,MUTED)
for x,color,desc in [(270,BLUE_FILL,'粗网格 h=2'),(635,FINE_FILL,'细网格 h=1'),(1000,OVER_FILL,'重叠计算区')]:
    layout.rect(x,152,28,24,color,GRID)
    layout.text(x+43,164,desc,28,anchor='start')
X0,Y0,S=270,225,960/1024
def X(x): return X0+x*S
def Y(y): return Y0+(ny-y)*S
def box(xlo,xhi,ylo,yhi,fill): layout.rect(X(xlo),Y(yhi),(xhi-xlo)*S,(yhi-ylo)*S,fill)
box(0,nx,0,ny,FINE_FILL)
box(wx,nx-wx,wy,ny-wy,BLUE_FILL)
# Decimated lines show the 2:1 grid-spacing relation; exact nodes are shown in the close-up.
for x in range(0,nx+1,32):
    for ya,yb in [(0,wy),(ny-wy,ny)]: layout.line(X(x),Y(ya),X(x),Y(yb),'#e9dcc8',1)
for y in range(0,ny+1,32):
    if y<=wy or y>=ny-wy: layout.line(X(0),Y(y),X(nx),Y(y),'#e9dcc8',1)
    else:
        layout.line(X(0),Y(y),X(wx),Y(y),'#e9dcc8',1)
        layout.line(X(nx-wx),Y(y),X(nx),Y(y),'#e9dcc8',1)
for x in list(range(0,wx+1,32))+list(range(nx-wx,nx+1,32)):
    layout.line(X(x),Y(wy),X(x),Y(ny-wy),'#e9dcc8',1)
for x in range(wx,nx-wx+1,64): layout.line(X(x),Y(wy),X(x),Y(ny-wy),'#cddbef',1)
for y in range(wy,ny-wy+1,64): layout.line(X(wx),Y(y),X(nx-wx),Y(y),'#cddbef',1)
# Shared overlap is twice the extension width at a straight interface.
for y in (wy,ny-wy): box(0,nx,y-ov+.5,y+ov+.5,OVER_FILL)
for x in (wx,nx-wx): box(x-ov+.5,x+ov+.5,wy-ov+.5,ny-wy+ov+.5,OVER_FILL)
# Physical domain and non-overlapping ownership boundaries.
layout.rect(X(0),Y(ny),nx*S,ny*S,None,INK,3)
for y in (wy,ny-wy): layout.line(X(0),Y(y),X(nx),Y(y),INK,2)
for x in (wx,nx-wx): layout.line(X(x),Y(wy),X(x),Y(ny-wy),INK,2)
# Artificial boundaries of all five extended computation rectangles.
for xa,xb,ya,yb in [(124.5,900.5,124.5,900.5),(0,1024,0,132.5),(0,1024,892.5,1024),(0,132.5,124.5,900.5),(892.5,1024,124.5,900.5)]:
    if xa>0: layout.line(X(xa),Y(ya),X(xa),Y(yb),OVER,1,True)
    if xb<nx: layout.line(X(xb),Y(ya),X(xb),Y(yb),OVER,1,True)
    if ya>0: layout.line(X(xa),Y(ya),X(xb),Y(ya),OVER,1,True)
    if yb<ny: layout.line(X(xa),Y(yb),X(xb),Y(yb),OVER,1,True)
layout.text(X(512),Y(960),'3  上细块   h=1   |   积分节点 1024 × 128',30)
layout.text(X(512),Y(64),'2  下细块   h=1   |   积分节点 1024 × 128',30)
layout.text(X(64),Y(585),'4\n左细块\nh=1',30)
layout.text(X(960),Y(585),'5\n右细块\nh=1',30)
layout.rect(X(265),Y(675),494*S,252*S,BLUE_FILL)
layout.text(X(512),Y(638),'1  中心粗块',44)
layout.text(X(512),Y(578),'h = 2   ·   Δt = 2',33)
layout.text(X(512),Y(522),'物理范围 768 × 768',32)
layout.text(X(512),Y(471),'积分节点 385 × 385',32)
layout.text(X(512),Y(298),'左右细块各为 128 × 768 个积分节点',28,MUTED)
layout.text(X(512),Y(252),'实线：物理积分分区    虚线：人工边界节点行 / 列',26,MUTED)
for tick in (0,128,896,1024):
    layout.line(X(tick),Y(0),X(tick),Y(0)+10,INK,2)
    layout.text(X(tick),Y(0)+34,str(tick),28)
    layout.text(X(0)-25,Y(tick),str(tick),28,anchor='end')
layout.text(X(512),Y(0)+79,'x  /  最细格距',29)
layout.text(195,190,'y',30)
layout.arrow(X(0),1300,X(wx),1300)
layout.text(328,1340,'左右加密厚度 = 128',27)
layout.arrow(X(wx),1300,X(nx-wx),1300)
layout.text(X(512),1340,'中心宽度 = 768',29)
layout.arrow(185,Y(0),185,Y(wy))
layout.text(165,Y(75),'上下加密厚度\n=128',26,anchor='end')
layout.text(750,1404,'概览线条仅示意格距比例；重合点见放大图。积分按分区裁剪面积，不重复计入面积。',26,MUTED)
layout.save('multiblock-layout')

detail=Canvas(1500,1220,'左细块与中心粗块的对齐格点；蓝框内橙点表示同坐标的粗细重合节点。')
detail.text(750,46,'当前接口：碰撞前接收，粗细节点重合',44)
detail.text(750,100,'粗节点是细节点的子集；蓝框内橙点 = 同一物理坐标、两块各自存储',28,MUTED)
ox,oy,scale=200,280,34
def xx(x): return ox+(x-112.5)*scale
def yy(y): return oy+(520.5-y)*scale
detail.rect(xx(112.5),yy(520.5),15.5*scale,12*scale,FINE_FILL)
detail.rect(xx(128),yy(520.5),16.5*scale,12*scale,BLUE_FILL)
detail.rect(xx(124.5),yy(520.5),8*scale,12*scale,OVER_FILL)
# Lines pass through lattice nodes, not cell faces. Both levels share the .5 phase.
for k in range(113,133): detail.line(xx(k+.5),yy(508.5),xx(k+.5),yy(520.5),'#e4cbb8',1)
for k in range(509,520): detail.line(xx(112.5),yy(k+.5),xx(132.5),yy(k+.5),'#e4cbb8',1)
for k in range(124,145,2): detail.line(xx(k+.5),yy(508.5),xx(k+.5),yy(520.5),'#a9bfdf',1)
for k in range(508,521,2): detail.line(xx(124.5),yy(k+.5),xx(144.5),yy(k+.5),'#a9bfdf',1)
detail.line(xx(128),yy(508.5),xx(128),yy(520.5),INK,2)
detail.line(xx(128.5),yy(508.5),xx(128.5),yy(520.5),OVER,2,True)
for k in range(112,133):
    for l in range(508,521): detail.dot(xx(k+.5),yy(l+.5),4,ORANGE)
for k in range(124,145,2):
    for l in range(508,521,2): detail.rect(xx(k+.5)-8,yy(l+.5)-8,16,16,None,BLUE,2)
detail.text(xx(124.5),244,'粗块首列 x=124.5',26)
detail.text(xx(132.5),244,'细块末列 x=132.5',26)
detail.arrow(xx(124.5),191,xx(132.5),191,OVER)
detail.text(xx(128.5),155,'重叠节点范围 [124.5,132.5]，跨度 8：9 排细节点、5 排粗节点',28,OVER)
for tick in (112.5,124.5,128.5,132.5,144.5): detail.text(xx(tick),735,str(tick),27)
for tick in (508.5,514.5,520.5): detail.text(175,yy(tick),str(tick),26,anchor='end')
detail.text(160,235,'y',28); detail.text(1340,735,'x',28)
detail.dot(205,795,4,ORANGE); detail.text(226,795,'细节点 h=1',27,anchor='start')
detail.rect(575,787,16,16,None,BLUE,2); detail.text(607,795,'粗节点 h=2',27,anchor='start')
detail.rect(965,787,16,16,None,BLUE,2); detail.dot(973,795,4,ORANGE)
detail.text(997,795,'粗细重合点',27,anchor='start')
detail.line(205,847,250,847,INK,2); detail.text(267,847,'实线 x=128：物理积分分区边界',26,anchor='start')
detail.line(205,892,250,892,OVER,2,True); detail.text(267,892,'虚线 x=128.5：共同节点列，两侧各延伸 4',26,anchor='start')
detail.text(750,945,'单侧延伸 4 = overlapCells × refineRatio = 2 × 2；Δx粗 / Δx细 = Δt粗 / Δt细 = 2',27)
detail.rect(140,987,1220,133,'#f5f7fa',GRID)
detail.text(170,1017,'interfaceSkin = 2：每次重建人工边界的两层节点',28,anchor='start')
detail.text(170,1068,'粗块首两列：124.5、126.5     |     细块末两列：131.5、132.5',27,anchor='start')
detail.text(750,1170,'重合点直接取值并做矩重标定；非重合细点仍做空间插值。物理壁面半格距保持不变。',26,MUTED)
detail.save('multiblock-interface')

clock=Canvas(1500,1170,'一次粗步：先准备碰撞前缓冲，粗步为2、细步为1，在t+2同步。')
clock.text(750,48,'读懂一次粗步：数据先到，节点再碰撞',43)
clock.text(750,104,'refineRatio = 2    |    以下时间以一个细步为单位',29,MUTED)
stages=[
    (170,'入口','两级都在 t；初始交换或上次同步已准备好缓冲',
     '交换量来自迁移后的 f、g，供接下来的碰撞使用', '#f1f4f8'),
    (345,'① 粗块预测','粗块：碰撞 → 迁移 → 原温度推进，到达 t+2',
     '保留 t−2、t、t+2 三个粗时间层；来源点排除两层人工边界',BLUE_FILL),
    (520,'② 细块第一步','准备 t 的接口 → 碰撞、迁移及温度推进 → t+1',
     '同级读取 t 快照；跨级读取粗块 t 快照，不能直接用 t+2',FINE_FILL),
    (695,'③ 细块第二步','准备 t+1 的接口 → 碰撞、迁移及温度推进 → t+2',
     '同级读取 t+1 快照；跨级由三个粗时间层插值到 t+1',FINE_FILL),
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

blocks=[]
for name,xa,xb,ya,yb,h in [('1 中心',128,896,128,896,2),('2 下',0,1024,0,128,1),
                          ('3 上',0,1024,896,1024,1),('4 左',0,128,128,896,1),('5 右',896,1024,128,896,1)]:
    ext=[max(0,xa-ov),min(nx,xb+ov),max(0,ya-ov),min(ny,yb+ov)]
    first=[ext[0]+.5,ext[2]+.5]
    last=[ext[1]-.5 if ext[1]==nx else ext[1]+.5,ext[3]-.5 if ext[3]==ny else ext[3]+.5]
    counts=[int((last[k]-first[k])/h)+1 for k in range(2)]
    owned_counts=[]
    for axis,(lo,hi) in enumerate([(xa,xb),(ya,yb)]):
        weights=[max(0,min(hi,first[axis]+i*h+h/2)-max(lo,first[axis]+i*h-h/2)) for i in range(counts[axis])]
        assert sum(weights)==hi-lo
        assert abs(sum(w*(first[axis]+i*h) for i,w in enumerate(weights))-(hi*hi-lo*lo)/2)<1e-9
        owned_counts.append(sum(w>0 for w in weights))
    blocks.append({'block':name,'owned_faces':[xa,xb,ya,yb],'h':h,
                   'first_computed_node':first,'last_computed_node':last,
                   'integration_nodes':owned_counts,'computed_nodes':counts})
manifest={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'nx':nx,'ny':ny,
          'refineRatio':ratio,'fineLayerCellsLeftRight':wx,'fineLayerCellsBottomTop':wy,'overlapCells':overlap_cells,
          'node_phase':0.5,'extension_from_reference_node':ov,'shared_node_span':2*ov,
          'interfaceSkin':skin,'blocks':blocks,
          'computed_node_total':sum(b['computed_nodes'][0]*b['computed_nodes'][1] for b in blocks),
          'integration_sample_total':sum(b['integration_nodes'][0]*b['integration_nodes'][1] for b in blocks)}
(HERE/'layout_parameters.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
