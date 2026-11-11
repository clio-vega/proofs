"""Q584 C1: is the PY hexagon (and trap) a composite of KZJ unit-triangle vertices?
Unit-triangle model, labels {0,1,X} (X='10').  Labels read from the Q569 rhombi split along the
short diagonal (derived in scratch notes):
  up   triangle (dir0,dir2,dir4) in cyclic class of (0,0,0),(1,1,1),(1,X,0)   [+ K-up (X,X,X)]
  down triangle (dir1,dir3,dir5) in cyclic class of (0,0,0),(1,1,1),(0,1,X)   [+ K-dn (X,X,X)]
Star(v) at v=(0,0): 6 triangles, 6 spokes, 6 rim edges. Sum over spoke labels."""
import itertools, sys
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import TILES, D, tile_edges, edge_key
X='X'
def cyc(t): return {t, t[1:]+t[:1], t[2:]+t[:2]}
UP=set().union(cyc((0,0,0)),cyc((1,1,1)),cyc((1,X,0)))
DN=set().union(cyc((0,0,0)),cyc((1,1,1)),cyc((0,1,X)))
def ek(p,q):
    d=D.index((q[0]-p[0],q[1]-p[1])); return edge_key(p,d)
def tri_edges(kind,a,b):
    if kind=='U': P=[(a,b),(a+1,b),(a,b+1)]; return [ek(P[0],P[1]),ek(P[1],P[2]),ek(P[2],P[0])]   # dir0,dir2,dir4
    P=[(a+1,b),(a+1,b+1),(a,b+1)]; return [ek(P[0],P[1]),ek(P[1],P[2]),ek(P[2],P[0])]            # dir1,dir3,dir5
def star(v):
    a,b=v
    return [('U',a,b),('D',a-1,b),('U',a-1,b),('D',a-1,b-1),('U',a,b-1),('D',a,b-1)]
def partition(tris, boundary, upK=0, dnK=0):
    """dict weightmonomial->count; boundary: edge->label.  K pieces weight tracked as exponents."""
    E=[e for t in tris for e in tri_edges(*t)]
    internal=sorted({e for e in E if E.count(e)==2})
    out={}
    for lab in itertools.product((0,1,X),repeat=len(internal)):
        L=dict(boundary); L.update(zip(internal,lab))
        ku=kd=0; ok=True
        for t in tris:
            trip=tuple(L[e] for e in tri_edges(*t))
            S=UP if t[0]=='U' else DN
            if trip in S: continue
            if trip==(X,X,X) and ((t[0]=='U' and upK) or (t[0]=='D' and dnK)):
                if t[0]=='U': ku+=1
                else: kd+=1
                continue
            ok=False; break
        if ok: out[(ku,kd)]=out.get((ku,kd),0)+1
    return out
if __name__=='__main__':
    # sanity / positive control: a KTW rhombus region (2 triangles) with rh_014 boundary must give 1
    anc=(0,0)
    B={e:l for e,l in tile_edges(TILES['rh_014'],anc)}
    print('control rh_014:',partition([('U',0,0),('D',0,0)],B))
    v=(0,0)
    # hex placed so its star centre is v: hex verts at anchor A are A+(0,0),(0,1),(-1,2),(-2,2),(-2,1),(-1,0); centre A+(-1,1)
    for nm in ('hex','hexR'):
        B={e:l for e,l in tile_edges(TILES[nm],(1,-1))}
        Es={e for t in star(v) for e in tri_edges(*t)}
        assert set(B)<=Es and len(B)==6, 'hex not on star'
        for uk,dk in ((0,0),(1,0),(0,1),(1,1)):
            print(nm,'K-up' if uk else '', 'K-dn' if dk else '', partition(star(v),B,uk,dk))
