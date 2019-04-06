import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

w = 96
h = 24

xdelta, ydelta = 32.7, 3
x,y = (
    np.arange(0, xdelta*4+xdelta, xdelta),
    -np.arange(0,ydelta*4+ydelta, ydelta)
)
boxes = []
for _x,_y in zip(x,y):
    rect = Rectangle((_x,_y), w, h, fill=True)
    boxes.append(rect)

pc = PatchCollection(boxes, facecolor='g', alpha=0.5, edgecolor='g',
                     linewidth=0)

f,ax = plt.subplots(figsize=(6,1.5))
ax.add_collection(pc)

ax.set_xlim((-5, np.max(x)+w+5))
ax.set_ylim((np.min(y)-5,np.max(y)+h+5))

plt.axis('off')
f.tight_layout(pad=0)

print('max x: {}'.format(np.max(x)))
print('max x + 96: {}'.format(np.max(x)+w))

outpath = '../results/rectangles_ecliptic.png'
f.savefig(outpath,dpi=350,bbox_inches='tight')
print('made {}'.format(outpath))
