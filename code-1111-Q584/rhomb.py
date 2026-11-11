"""Q584 C3: charge-conserving rhombi; can hex/hexR rim be filled by 3 of them (lozenge tilings of star),
or by mixed rhombus+triangle fillings?  Support argument: every U_q(A2) intertwiner entry conserves weight."""
from charge import N, idx, NC, row
from star import tri_edges, star, UP, DN, X, ek
import itertools, sys
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import TILES, tile_edges, D, edge_key
Nn=[v for v in N if v[NC]==0]          # charges with c=0 (rhombi have 2 tris of opposite orientation -> c cancels)
def conserved(edges):   # edges: list of (dir,label) ccw
    r=row(edges,0); return all(sum(a*b for a,b in zip(r,v))==0 for v in N)
# rhombus shapes: ccw step dirs
SH={'r01':(0,1,3,4),'r12':(1,2,4,5),'r20':(2,3,5,0)}
RH={}
for nm,dirs in SH.items():
    RH[nm]=[lab for lab in itertools.product((0,1,X),repeat=4) if conserved(list(zip(dirs,lab)))]
    print(nm,'conserving labelings:',len(RH[nm]))
# positive control: KTW rhombi must be conserving
for t in ('rh_014','rh_230','rh_125'): print('control',t,conserved(TILES[t]))
# generic region filler: pieces = unit triangles (UP/DN + XXX) and conserving rhombi, exact cover of star(v)
v=(0,0); TRI=star(v)
def tri_key(t): return t
def rh_place():
    # all rhombi = pairs of adjacent triangles in star sharing an edge
    out=[]
    for s,t in itertools.combinations(TRI,2):
        sh=set(tri_edges(*s))&set(tri_edges(*t))
        if sh:
            e=sh.pop()
            outer=[x for x in tri_edges(*s)+tri_edges(*t) if x!=e]
            out.append(((s,t),outer))
    return out
RP=rh_place()
def dir_of(e,cc_poly):
    return None
def fill(rim, allow_rh=True, allow_tri=True):
    """count fillings; rhombus validity via charge conservation on its 4 boundary edges, using
    orientation: we evaluate conservation via edge keys (class, sign) relative to the piece interior."""
    res=0
    # brute: assign labels to 6 spokes, then check all decompositions into pieces
    rimL=dict(rim)
    spokes=sorted({e for t in TRI for e in tri_edges(*t)}-set(rimL))
    W=[]
    for lab in itertools.product((0,1,X),repeat=len(spokes)):
        L=dict(rimL); L.update(zip(spokes,lab))
        # count exact covers of TRI by pieces valid under L
        def ok_tri(t):
            trip=tuple(L[e] for e in tri_edges(*t)); return trip in (UP if t[0]=='U' else DN) or trip==(X,X,X)
        def ok_rh(pair,outer):
            # build ccw (dir,label) list of the rhombus boundary: orient each edge with interior on left
            s,t=pair; es=[]
            for tt in (s,t):
                kind,a,b=tt
                P=[(a,b),(a+1,b),(a,b+1)] if kind=='U' else [(a+1,b),(a+1,b+1),(a,b+1)]
                for i in range(3):
                    p,q=P[i],P[(i+1)%3]; d=D.index((q[0]-p[0],q[1]-p[1]))
                    if edge_key(p,d) in outer: es.append((d,L[edge_key(p,d)]))
            return conserved(es)
        def covers(rem):
            if not rem: return 1
            t=min(rem); c=0
            if allow_tri and ok_tri(t): c+=covers(rem-{t})
            if allow_rh:
                for pair,outer in RP:
                    if t in pair and set(pair)<=rem and ok_rh(pair,outer):
                        c+=covers(rem-set(pair))
            return c
        c=covers(frozenset(TRI))
        if c: W.append((lab,c))
    return W
rim={e for t in TRI for e in tri_edges(*t) if sum(e in tri_edges(*s) for s in TRI)==1}
def place(nm):
    for a in range(-3,4):
        for b in range(-3,4):
            B={e:l for e,l in tile_edges(TILES[nm],(a,b))}
            if set(B)==rim: return B
for nm in ('hex','hexR'):
    B=place(nm)
    print(nm,'triangles only   :',len(fill(B,False,True)))
    w=fill(B,True,False); print(nm,'rhombi only      :',len(w), w[:4])
    w=fill(B,True,True);  print(nm,'mixed            :',len(w))
