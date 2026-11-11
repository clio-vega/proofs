from star import *
v=(0,0); Es={e for t in star(v) for e in tri_edges(*t)}
rim={e for t in star(v) for e in tri_edges(*t) if sum(e in tri_edges(*s) for s in star(v))==1}
def place(nm):
    for a in range(-3,4):
        for b in range(-3,4):
            B={e:l for e,l in tile_edges(TILES[nm],(a,b))}
            if set(B)==rim: return B
for nm in ('hex','hexR'):
    B=place(nm); print(nm,'rim labels found:',B is not None)
    for uk,dk in ((0,0),(1,0),(0,1),(1,1)):
        print('  ',nm,'upK' if uk else '-', 'dnK' if dk else '-', partition(star(v),B,uk,dk))
# what rim labellings in {0,1}^6 DO admit fillings, and with which K usage?
import itertools
rimL=sorted(rim)
res={}
for lab in itertools.product((0,1),repeat=6):
    p=partition(star(v),dict(zip(rimL,lab)),1,1)
    if p: res[lab]=p
print(len(res),'of 64 {0,1}-rims fillable');
for k,p in res.items():
    if any(m!=(0,0) for m in p): print('  K-using rim',k,p)
