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
nx = int(re.search(r'^#define NX_OVERRIDE (\d+)', src, re.M)[1])
ny = int(re.search(r'^#define NY_OVERRIDE (\d+)', src, re.M)[1])
ratio = int(re.search(r'parameter :: refineRatio=(\d+)', src)[1])
overlap_cells = int(re.search(r'parameter :: overlapCells=(\d+)', src)[1])
skin = int(re.search(r'parameter :: interfaceSkin=(\d+)', src)[1])
assert 'wallCellsX=nx/8, wallCellsY=ny/8' in src
wx, wy, ov = nx // 8, ny // 8, overlap_cells * ratio
assert (nx, ny, ratio, wx, wy, ov, skin) == (1024, 1024, 2, 128, 128, 8, 2), \
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
for y in (wy,ny-wy): box(0,nx,y-ov,y+ov,OVER_FILL)
for x in (wx,nx-wx): box(x-ov,x+ov,wy-ov,ny-wy+ov,OVER_FILL)
# Physical domain and non-overlapping ownership boundaries.
layout.rect(X(0),Y(ny),nx*S,ny*S,None,INK,3)
for y in (wy,ny-wy): layout.line(X(0),Y(y),X(nx),Y(y),INK,2)
for x in (wx,nx-wx): layout.line(X(x),Y(wy),X(x),Y(ny-wy),INK,2)
# Artificial boundaries of all five extended computation rectangles.
for xa,xb,ya,yb in [(120,904,120,904),(0,1024,0,136),(0,1024,888,1024),(0,136,120,904),(888,1024,120,904)]:
    if xa>0: layout.line(X(xa),Y(ya),X(xa),Y(yb),OVER,1,True)
    if xb<nx: layout.line(X(xb),Y(ya),X(xb),Y(yb),OVER,1,True)
    if ya>0: layout.line(X(xa),Y(ya),X(xb),Y(ya),OVER,1,True)
    if yb<ny: layout.line(X(xa),Y(yb),X(xb),Y(yb),OVER,1,True)
layout.text(X(512),Y(960),'3  上细块   h=1   |   所有权节点 1024 × 128',30)
layout.text(X(512),Y(64),'2  下细块   h=1   |   所有权节点 1024 × 128',30)
layout.text(X(64),Y(585),'4\n左细块\nh=1',30)
layout.text(X(960),Y(585),'5\n右细块\nh=1',30)
layout.rect(X(265),Y(675),494*S,252*S,BLUE_FILL)
layout.text(X(512),Y(638),'1  中心粗块',44)
layout.text(X(512),Y(578),'h = 2   ·   Δt = 2',33)
layout.text(X(512),Y(522),'物理范围 768 × 768',32)
layout.text(X(512),Y(471),'所有权节点 384 × 384',32)
layout.text(X(512),Y(298),'左右细块各为 128 × 768 个所有权节点',28,MUTED)
layout.text(X(512),Y(252),'实线：不重复计数的分区边界    虚线：计算区边界',27,MUTED)
for tick in (0,128,896,1024):
    layout.line(X(tick),Y(0),X(tick),Y(0)+10,INK,2)
    layout.text(X(tick),Y(0)+34,str(tick),28)
    layout.text(X(0)-25,Y(tick),str(tick),28,anchor='end')
layout.text(X(512),Y(0)+79,'x  /  最细格距',29)
layout.text(195,190,'y',30)
layout.arrow(X(0),1300,X(wx),1300)
layout.text(328,1340,'wallCellsX = 128',27)
layout.arrow(X(wx),1300,X(nx-wx),1300)
layout.text(X(512),1340,'nx − 2×wallCellsX = 768',29)
layout.arrow(185,Y(0),185,Y(wy))
layout.text(165,Y(75),'wallCellsY\n=128',26,anchor='end')
layout.text(750,1404,'概览网格线抽稀 32 倍；色带为实际比例的重叠区，统计与输出只计各块所有权区。',27,MUTED)
layout.save('multiblock-layout')

detail=Canvas(1500,1170,'左细块与中心粗块的真实格点；每块延伸8，共同覆盖16，人工边界分别重建两层节点。')
detail.text(750,46,'左细块与中心粗块：接口放大',44)
detail.text(750,102,'真实节点布置    ·    分区接口 x = 128    ·    所有长度以最细格距计',28,MUTED)
ox,oy,scale=260,250,30
def xx(x): return ox+(x-112)*scale
def yy(y): return oy+(520-y)*scale
detail.rect(xx(112),yy(520),16*scale,16*scale,FINE_FILL)
detail.rect(xx(128),yy(520),16*scale,16*scale,BLUE_FILL)
detail.rect(xx(120),yy(520),16*scale,16*scale,OVER_FILL)
for x in range(112,137): detail.line(xx(x),yy(504),xx(x),yy(520),'#ead5be',1)
for y in range(504,521): detail.line(xx(112),yy(y),xx(136),yy(y),'#ead5be',1)
for x in range(120,145,2): detail.line(xx(x),yy(504),xx(x),yy(520),'#b7cbea',2)
for y in range(504,521,2): detail.line(xx(120),yy(y),xx(144),yy(y),'#b7cbea',2)
for x in range(112,136):
    for y in range(504,520): detail.dot(xx(x+.5),yy(y+.5),3.7,ORANGE)
for x in range(121,144,2):
    for y in range(505,520,2): detail.dot(xx(x),yy(y),5,BLUE,True)
detail.rect(xx(112),yy(520),32*scale,16*scale,None,GRID,2)
detail.line(xx(128),yy(504),xx(128),yy(520),INK,3)
detail.line(xx(120),yy(504),xx(120),yy(520),BLUE,3,True)
detail.line(xx(136),yy(504),xx(136),yy(520),ORANGE,3,True)
detail.text(xx(120),225,'粗块计算区起点',27)
detail.text(xx(136),225,'细块计算区终点',27)
detail.arrow(xx(120),177,xx(136),177,OVER)
detail.text(xx(128),148,'共同覆盖宽度 = 16',29,OVER)
for tick in (112,120,128,136,144): detail.text(xx(tick),762,str(tick),28)
for tick in (504,512,520): detail.text(235,yy(tick),str(tick),27,anchor='end')
detail.text(235,215,'y',28)
detail.text(1260,762,'x',28)
detail.dot(310,795,4,ORANGE); detail.text(328,795,'细节点 h=1',26,anchor='start')
detail.dot(1050,795,5,BLUE,True); detail.text(1068,795,'粗节点 h=2',26,anchor='start')
detail.arrow(xx(120),815,xx(128),815,BLUE)
detail.arrow(xx(128),815,xx(136),815,ORANGE)
detail.text(xx(124),851,'粗块向左延伸 8',26)
detail.text(xx(132),851,'细块向右延伸 8',26)
detail.text(750,890,'8 = overlapCells × refineRatio = 4 × 2',26)
detail.text(750,929,'interfaceSkin = 2：人工边界重建的两层节点',29)
detail.text(235,986,'粗块覆盖层',27,anchor='end')
detail.rect(xx(120),972,24*scale,28,BLUE_FILL)
detail.rect(xx(120),972,4*scale,28,BLUE)
detail.text(650,986,'x∈[120,124]，2 个粗格，宽度 4',26,anchor='start')
detail.text(235,1050,'细块覆盖层',27,anchor='end')
detail.rect(xx(112),1036,24*scale,28,FINE_FILL)
detail.rect(xx(134),1036,2*scale,28,ORANGE)
detail.text(300,1050,'x∈[134,136]，2 个细格，宽度 2',26,anchor='start')
detail.text(750,1125,'overlapCells 控制计算区延伸；interfaceSkin 控制人工边界覆盖，两者含义不同。',27,MUTED)
detail.save('multiblock-interface')

blocks=[]
for name,xa,xb,ya,yb,h in [('1 中心',128,896,128,896,2),('2 下',0,1024,0,128,1),
                          ('3 上',0,1024,896,1024,1),('4 左',0,128,128,896,1),('5 右',896,1024,128,896,1)]:
    comp=[max(0,xa-ov),min(nx,xb+ov),max(0,ya-ov),min(ny,yb+ov)]
    blocks.append({'block':name,'owned_faces':[xa,xb,ya,yb], 'computed_faces':comp,'h':h,
                   'owned_nodes':[(xb-xa)//h,(yb-ya)//h],
                   'computed_nodes':[(comp[1]-comp[0])//h,(comp[3]-comp[2])//h]})
manifest={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'nx':nx,'ny':ny,
          'refineRatio':ratio,'wallCellsX':wx,'wallCellsY':wy,'overlapCells':overlap_cells,
          'extension_fine_units':ov,'shared_overlap_fine_units':2*ov,'interfaceSkin':skin,'blocks':blocks}
(HERE/'layout_parameters.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
