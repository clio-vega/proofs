"""PROVE 11-21: the two-word model for (CS) -- is the first contact of y's sweep with Jo's net-removed set W always
y's straight step (class 1; a1 for y->B, a0 for y->C)?
Model: Jo (target T_J) from the nest-A rhombus at nest position x, y (target T_y != T_J) from nest position y, both
FREE words in {a0,a1,b} (each step = the unique good configuration of its (class,type), 11-20 Lemma word).
Constraints (necessary for both journeys to run on one start board M0):
  (K1) every tile a journey removes that it did not itself create earlier is a tile of M0; all such tiles of both
       journeys, together with the nest-A rhombi of the minimal ideal containing x and y, are pairwise EQUAL or
       interior-disjoint (plane polygons, SAT test with tolerance);
  (K2) every such tile lies in the 150-degree nest cone at A (Purbhoo mosaic, Lemma N).
W := M0-tiles removed by Jo and absent at the end of Jo's word.  Contact := first y-step removing a tile of W.
The engine's choice rule is NOT modelled, so a violation here may be unrealisable; zero violations is evidence.
twoword2 = twoword + contact pairing (y class,type) vs (Jo class,type) at removal.  usage: python3 twoword2.py LJ LY [margin]"""
import sys, cmath, itertools, os
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597')
sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597')
from fractions import Fraction as Fr
from collections import Counter
import words as Wd
import localenum as LE
LJ, LY = int(sys.argv[1]), int(sys.argv[2]); MARG = int(sys.argv[3]) if len(sys.argv) > 3 else 0
FIRSTA0 = os.environ.get('FIRSTA0', '0') == '1'   # constraint: both journeys start with a0 (4779/4779 on the dumps)
UA = (0, -1, 0, 0); WA = (0, 0, 0, 1)
def it(t): return (t[0], tuple(sorted(tuple(int(x) for x in P) for P in t[1])))
def tr(t, v): return (t[0], tuple(sorted(tuple(a+b for a, b in zip(P, v)) for P in t[1])))
CF = {}
for T in 'BC':
    C = Wd.configs(T)
    for (ck, typ), (Rc, rem, Rn, cre) in C.items():
        CF[(T, ck, typ)] = (it(Rc), [it(t) for t in rem], it(Rn), [it(t) for t in cre])
def classkey(R):
    d = set()
    for A in R[1]:
        for B in R[1]:
            w = tuple(b-a for a, b in zip(A, B))
            if any(w[:2]) and not any(w[2:]): d.add(w[:2])
    return frozenset(d)
def base(R): return min(R[1])
def step(T, R, typ):
    ck = classkey(R); Rc, rem, Rn, cre = CF[(T, ck, typ)]
    sh = base(R); assert tr(Rc, sh) == R
    return [tr(t, sh) for t in rem], tr(Rn, sh), [tr(t, sh) for t in cre]
# plane geometry
POLY = {}
def poly(t):
    if t not in POLY:
        zs = [LE.phi(P) for P in t[1]]; POLY[t] = LE.convex_hull(zs)
    return POLY[t]
OV = {}
def overlap(s, t):
    """True iff interiors intersect (SAT)."""
    key = (s, t) if s < t else (t, s)
    if key in OV: return OV[key]
    A, B = poly(s), poly(t); res = True
    for P in (A, B):
        for i in range(len(P)):
            e = P[(i+1) % len(P)] - P[i]; n = complex(-e.imag, e.real)
            pa = [(z*n.conjugate()).real for z in A]; pb = [(z*n.conjugate()).real for z in B]
            if max(pa) <= min(pb) + 1e-9 or max(pb) <= min(pa) + 1e-9: res = False; break
        if not res: break
    OV[key] = res; return res
zu, zw = LE.phi(UA), LE.phi(WA)
def incone(t):
    for P in t[1]:
        z = LE.phi(P)
        # z = s zu + t zw
        det = (zu.conjugate()*zw).imag
        s = (z.conjugate()*zw).imag / det * -1 if False else ((z*zw.conjugate()).imag) / ((zu*zw.conjugate()).imag)
        tt = ((z*zu.conjugate()).imag) / ((zw*zu.conjugate()).imag)
        if s < -1e-9 or tt < -1e-9: return False
    return True
def nestrh(i, j):
    p = tuple(i*a + j*b for a, b in zip(UA, WA))
    return ('rh', tuple(sorted(tuple(a+b for a, b in zip(p, d)) for d in [(0,)*4, UA, WA, tuple(a+b for a, b in zip(UA, WA))])))
def compatible(t, S):
    return all(t == s or not overlap(t, s) for s in S)
def strt(T): return {'B': (1, 'a1'), 'C': (1, 'a0')}[T]
st = Counter()
OFFS = {('B', 'C'): [(1, -1), (1, -2), (2, -1), (2, -2), (3, -1)], ('C', 'B'): [(-1, 1), (-2, 1), (-1, 2), (-2, 2), (-1, 3)]}
TYPES = ('a0', 'a1', 'b')
for (Ty, TJ), offs in OFFS.items():
    for (di, dj) in offs:
        # minimal ideal containing x and y, plus a margin
        ix = (max(0, -di) if True else 0) + MARG; jx = max(0, -dj) + MARG
        iy, jy = ix + di, jx + dj
        ideal = {(i, j) for i in range(ix+1) for j in range(jx+1)} | {(i, j) for i in range(iy+1) for j in range(jy+1)}
        NEST = [nestrh(i, j) for (i, j) in ideal]
        Rx, Ry = nestrh(ix, jx), nestrh(iy, jy)
        nestother = [r for r in NEST if r not in (Rx, Ry)]
        # all Jo words, with consistency
        REMAT = {}
        def jdfs(R, read, created, word):
            if word:
                fin = set(created)
                W = frozenset(t for t in read if t not in fin)
                yield word, read, W, dict(REMAT)
            if len(word) == LJ: return
            for typ in (('a0',) if (FIRSTA0 and not word) else TYPES):
                rem, Rn, cre = step(TJ, R, typ)
                newread = list(read); ok = True
                for t in rem:
                    if t[0] == 'rh': continue
                    if t in created: continue
                    if not incone(t) or not compatible(t, nestother + newread): ok = False; break
                    newread.append(t)
                if not ok: st[(Ty, 'Jo pruned')] += 1; continue
                nc = (set(created) - set(rem)) | set(c for c in cre if c[0] != 'rh')
                cJ = 0 if classkey(R) == classkey(Rx) else 1
                added = [t for t in rem if t[0] != 'rh' and t not in REMAT]
                for t in added: REMAT[t] = (cJ, typ)
                yield from jdfs(Rn, newread, nc, word + [typ])
                for t in added: del REMAT[t]
        for jw, jread, W, RA in jdfs(Rx, [], set(), []):
            if not W: continue
            st[(Ty, 'Jo words')] += 1
            base_ctx = nestother + list(jread)
            def ydfs(R, read, created, word, cls):
                if len(word) == LY: st[(Ty, 'no contact within LY')] += 1; return
                for typ in (('a0',) if (FIRSTA0 and not word) else TYPES):
                    rem, Rn, cre = step(Ty, R, typ)
                    ck0 = cls
                    hit = [t for t in rem if t[0] != 'rh' and t not in created and t in W]
                    newread = list(read); ok = True
                    for t in rem:
                        if t[0] == 'rh' or t in created: continue
                        if t in W: continue
                        if not incone(t) or not compatible(t, base_ctx + newread): ok = False; break
                        newread.append(t)
                    # Jo's start rhombus is a tile of M0 too
                    if ok and any(t == Rx for t in rem): ok = False
                    if not ok: continue
                    c = 0 if classkey(R) == classkey(Ry) else 1
                    if hit:
                        good = (c, typ) == strt(Ty)
                        st[(Ty, 'CONTACT', 'straight' if good else 'NOT straight', (c, typ), 'kinds', tuple(sorted(t[0] for t in hit)))] += 1
                        st[(Ty, 'PAIR y', (c, typ), 'Jo at removal', tuple(sorted((t[0],) + RA[t] for t in hit)))] += 1
                        if not good and st[(Ty, 'shown')] < 3:
                            st[(Ty, 'shown')] += 1
                            print('VIOL', Ty, (di, dj), 'Jo', ''.join(x[-1] if x != 'b' else 'b' for x in jw), jw, 'y', word + [typ], 'class', c, flush=True)
                        continue
                    nc = (set(created) - set(rem)) | set(t for t in cre if t[0] != 'rh')
                    ydfs(Rn, newread, nc, word + [typ], cls)
            ydfs(Ry, [], set(), [], 0)
    print('done', Ty, flush=True)
for k in sorted(st, key=str): print(st[k], k)
