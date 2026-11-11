"""Independent check (edge-set, not traversal): sigma_W = Euclidean reflection in a=0,
linear map (a,b)->(-a,a+b) [fixes e1, sends e0->e2]. Compare labelled edge SETS."""
import sys; sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import TILES, D, tile_edges, verts
def sig(p): return (-p[0], p[0]+p[1])
def edge_pts(ek):
    p,d=ek; q=(p[0]+D[d][0],p[1]+D[d][1]); return frozenset([p,q])
def labelled(nm,anchor): return {edge_pts(e):l for e,l in tile_edges(TILES[nm],anchor)}
def centre(nm,anchor):
    V=verts(TILES[nm],anchor); return (sum(v[0] for v in V), sum(v[1] for v in V), len(V))
C={0:1,1:0,2:2}
anc=(0,0)
H=labelled('hex',anc)
img={frozenset(sig(p) for p in e):l for e,l in H.items()}
imgC={e:C[l] for e,l in img.items()}
# find anchors of tiles whose labelled edge set equals img / imgC
for nm in ('hex','hexR'):
  for a in range(-4,5):
    for b in range(-4,5):
      L=labelled(nm,(a,b))
      if L==img: print('sigma_W(hex) =',nm,'at',(a,b))
      if L==imgC: print('sigma_W∘C(hex) =',nm,'at',(a,b))
# centre of hex at anchor (0,0) vs its image: is the hex centred on the wall?
print('hex vertices',verts(TILES['hex'],anc))
# negative/positive control: identity map must return hex itself
print('control id:', any(labelled('hex',(a,b))==H for a in range(-2,3) for b in range(-2,3)))
