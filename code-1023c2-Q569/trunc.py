from puzzles import *

def raw_cover(steps, anchor):
    poly=verts(steps,anchor); area2=ccw_area2(steps); got=[]
    for a in range(-4,6):
        for b in range(-4,8):
            for tag,tri in (("U",up_tri(a,b)),("D",dn_tri(a,b))):
                if inside(poly,tri_centroid3(tri)): got.append((tag,a,b))
    assert len(got)==area2,(len(got),area2)
    return set(got)

def tri_minA(tag,a,b):
    t = up_tri(a,b) if tag=="U" else dn_tri(a,b)
    return min(v[0] for v in t)

b0=0
H  = raw_cover(TILES['hex'], (1,b0))     # hexagon anchored one step too far left
T  = raw_cover(TILES['trap'],(1,b0))
print("hexagon anchored at (1,0): covers", sorted(H))
print("  centroid of hexagon (a-coord, x3) =", sum(v[0] for v in verts(TILES['hex'],(1,b0))))
print("  -> a-coords of its 6 vertices:", [v[0] for v in verts(TILES['hex'],(1,b0))])
Hin  = {t for t in H if tri_minA(*t)>=0}
Hout = {t for t in H if tri_minA(*t)< 0}
print()
print("hexagon  INSIDE the wall {a>=0} :", sorted(Hin))
print("hexagon OUTSIDE the wall        :", sorted(Hout))
print("trapezoid covers                :", sorted(T))
print()
print("  ** trap == hex INTERSECT {a>=0} :", T==Hin, "**")
print("  halves equal in size           :", len(Hin)==len(Hout)==3)
eH={e:l for e,l in tile_edges(TILES['hex'],(1,b0))}
eT={e:l for e,l in tile_edges(TILES['trap'],(1,b0))}
shared={e:l for e,l in eH.items() if e in eT}
print()
print("  hex edges kept by trap, with labels:", {k:(eH[k],eT[k]) for k in shared}, "-> labels agree:", all(eH[k]==eT[k] for k in shared))
print("  NEW edges of trap (the cut, on the wall a=0):", {e:l for e,l in eT.items() if e not in eH})
print("  hex edges DISCARDED:", {e:l for e,l in eH.items() if e not in eT})
