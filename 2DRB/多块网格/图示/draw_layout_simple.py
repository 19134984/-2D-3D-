"""从当前源码读取分区比例，输出只含粗细网格标注的示意图。"""
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
source = (HERE.parent / '2DRBOpenaccMultiblock.F90').read_text(encoding='utf-8-sig')
source = '\n'.join(line.split('!')[0] for line in source.splitlines())

def parameter(name):
    return int(re.search(r'\b' + name + r'\s*=\s*(\d+)', source, re.I)[1])

nx, ny, ratio = map(parameter, ('nx', 'ny', 'refineRatio'))
left, right, bottom, top = [parameter('fineLayerCells' + side)
                            for side in ('Left', 'Right', 'Bottom', 'Top')]
W = H = 1200
S = 2
im = Image.new('RGB', (W*S, H*S), 'white')
draw = ImageDraw.Draw(im)
svg = ET.Element('svg', xmlns='http://www.w3.org/2000/svg', width=str(W),
                 height=str(H), viewBox=f'0 0 {W} {H}')

def rect(x, y, w, h, fill, stroke=None, width=1):
    draw.rectangle(tuple(round(v*S) for v in (x,y,x+w,y+h)),
                   fill=fill, outline=stroke, width=width*S)
    ET.SubElement(svg, 'rect', x=str(x), y=str(y), width=str(w), height=str(h),
                  fill=fill or 'none', stroke=stroke or 'none',
                  **{'stroke-width': str(width)})

def line(x1, y1, x2, y2, color, width=1):
    draw.line(tuple(round(v*S) for v in (x1,y1,x2,y2)), fill=color, width=width*S)
    ET.SubElement(svg, 'line', x1=str(x1), y1=str(y1), x2=str(x2), y2=str(y2),
                  stroke=color, **{'stroke-width': str(width)})

def label(x, y, text, color, fill):
    rect(x-113, y-36, 226, 72, fill)
    font = ImageFont.truetype('C:/Windows/Fonts/msyh.ttc', 44*S)
    draw.text((round(x*S), round(y*S)), text, font=font, fill=color, anchor='mm')
    el = ET.SubElement(svg, 'text', x=str(x), y=str(y), fill=color,
                       **{'font-family': 'Microsoft YaHei, sans-serif',
                          'font-size': '44', 'text-anchor': 'middle',
                          'dominant-baseline': 'central'})
    el.text = text

rect(0, 0, W, H, '#ffffff')
x0 = y0 = 80
extent = 1040
x1 = x0 + extent
y1 = y0 + extent
xl = x0 + extent*(left-.5)/nx
xr = x0 + extent*(nx-right+.5)/nx
yt = y0 + extent*(top-.5)/ny
yb = y0 + extent*(ny-bottom+.5)/ny
fine_fill, coarse_fill = '#fff4e6', '#edf3fc'
rect(x0, y0, extent, extent, fine_fill)
rect(xl, yt, xr-xl, yb-yt, coarse_fill)

# 适量抽样显示网格线，保留实际粗细格距比，避免全量节点挤成色块。
step = extent/48
for i in range(49):
    x = x0 + i*step
    if xl < x < xr:
        line(x,y0,x,yt,'#d9b98d')
        line(x,yb,x,y1,'#d9b98d')
    else:
        line(x,y0,x,y1,'#d9b98d')
    y = y0 + i*step
    if yt < y < yb:
        line(x0,y,xl,y,'#d9b98d')
        line(xr,y,x1,y,'#d9b98d')
    else:
        line(x0,y,x1,y,'#d9b98d')

coarse_step = step*ratio
i = 1
while xl+i*coarse_step < xr:
    line(xl+i*coarse_step,yt,xl+i*coarse_step,yb,'#91aacb')
    i += 1
i = 1
while yt+i*coarse_step < yb:
    line(xl,yt+i*coarse_step,xr,yt+i*coarse_step,'#91aacb')
    i += 1
rect(xl,yt,xr-xl,yb-yt,None,'#5179ac',3)
rect(x0,y0,extent,extent,None,'#35465b',3)
label((xl+xr)/2, (yt+yb)/2, '粗网格', '#315b8b', coarse_fill)
label((x0+x1)/2, (y0+yt)/2, '细网格', '#9a5f20', fine_fill)

im.save(HERE / 'multiblock-layout-simple.png', dpi=(300,300))
ET.indent(svg)
ET.ElementTree(svg).write(HERE / 'multiblock-layout-simple.svg',
                          encoding='utf-8', xml_declaration=True)
assert [el.text for el in svg.iter('text')] == ['粗网格', '细网格']
print('Saved multiblock-layout-simple.png and multiblock-layout-simple.svg')
