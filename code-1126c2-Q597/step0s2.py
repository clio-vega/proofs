"""PROVE 2026-11-26 c2, Step 0: (Share_2) as a predicate, its translation classes, and DREAM's zone falsifier.
Reads the record pickles written by code-1126-Q597/step0.py: /tmp/step0_{SIZE}_*.pkl.
Record = (case, L, [(T, A, Ar, p, |midA|, |midAr|) for the two roles]).  Predicates: see PREDICTION.md.
usage: python3 step0s2.py SIZE"""
import sys, glob, pickle
from collections import Counter
SIZE = int(sys.argv[1])
UNITS = {tuple(int(i == j) * sg for j in range(4)) for i in range(4) for sg in (1, -1)} | \
        {tuple(sg*c for c in u) for u in ((1, -1, 0, 0), (0, 0, 1, -1)) for sg in (1, -1)}
def vsub(P, Q): return tuple(a - b for a, b in zip(P, Q))
def G(s): return {t for t in s[1] if t[0] != 'rh'}
def rhomb(S): return next(t for t in S if t[0] == 'rh')
def trT(t, v): return (t[0], frozenset(tuple(a - b for a, b in zip(P, v)) for P in t[1]))
def canon(s):
    """(shape, base point) of a step: translation class and offset"""
    pts = [P for t in list(s[1]) + list(s[2]) for P in t[1]]; b = min(pts)
    return ((s[0], frozenset(trT(t, b) for t in s[1]), frozenset(trT(t, b) for t in s[2])), b)
def edges(t):
    P = list(t[1]); return {frozenset((A, B)) for A in P for B in P if vsub(B, A) in UNITS}
def alpha_edges(t, al): return {e for e in edges(t) if vsub(*sorted(e, key=lambda X: X)) in (al, tuple(-x for x in al))}
st = Counter(); classes = Counter(); ex = {}
for fn in sorted(glob.glob('/tmp/step0_%d_*.pkl' % SIZE)):
    for case, L, info in pickle.load(open(fn, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        E0, E1, L0, L1 = AE[pE], AE[pE+1], AL[pL], AL[pL+1]
        HE, HL = set(E0[1]), set(L0[1]); CE0, CL0 = set(E0[2]), set(L0[2])
        u = [t for t in HE & HL if t[0] == 'sq']
        st[(case, 'shared square exists (Thm share)', len(u) == 1)] += 1
        if len(u) != 1: continue
        u = u[0]
        GE1, GL1 = G(E1), G(L1)
        # (S2.i)
        iE = (not (GE1 & CE0)) and len(GE1 & HL) == 1 and not (GE1 & HE)
        iL = (not (GL1 & CL0)) and len(GL1 & HE) == 1 and not (GL1 & HL)
        st[(case, 'S2.i  E side', iE, 'L side', iL)] += 1
        st[(case, 'E1 rhombus created by E0', rhomb(E1[1]) in CE0, 'L1 rhombus by L0', rhomb(L1[1]) in CL0)] += 1
        OE, OL = GE1 - HE - HL, GL1 - HE - HL
        st[(case, 'S2.ii Oth(E1)==Oth(L1)', OE == OL, '|Oth|', len(OE), len(OL))] += 1
        st[(case, 'Oth kinds', tuple(sorted(t[0] for t in OE)), 'hexL tile of E1', tuple(t[0] for t in GE1 & HL),
            'hexE tile of L1', tuple(t[0] for t in GL1 & HE))] += 1
        cc = set(E1[2]) & set(L1[2])
        st[(case, 'S2.iii |C(E1)&C(L1)|', len(cc), 'kinds', tuple(sorted(t[0] for t in cc)))] += 1
        # translation class of the joint tuple, relative to the shared square
        b = min(u[1])
        key = (case, TL, tuple(frozenset(trT(t, b) for t in S) for S in (E0[1], E0[2], E1[1], E1[2], L0[1], L0[2], L1[1], L1[2])))
        classes[key] += 1
        ex.setdefault(key, (E0[0], E1[0], L0[0], L1[0]))
        # zone falsifier, case b with corridor (k>0)
        k = laL - 3 if case == 'b' else 0
        if case == 'b' and k > 0:
            (a1, m1), (b1, n1) = canon(AL[pL+1]), canon(ArL[pL])
            if a1 != b1: st[('zone', 'L^E_0 not a translate of L_1')] += 1; continue
            al = vsub(m1, n1)   # L^E_0 = L_1 - alpha
            st[('zone', 'alpha is a unit vector', al in UNITS)] += 1
            for name, cor in (('E^L', [ArE[pE+i] for i in range(2, k+2)]), ('L^E', [ArL[pL+i] for i in range(1, k+1)]),
                              ('J_E (M0)', [AE[pE+i] for i in range(2, k+2)])):
                tiles = [t for s in cor for t in s[1]]
                has = [bool(alpha_edges(t, al)) for t in tiles]
                st[('zone', name, 'every corridor tile has an alpha-edge', all(has), 'frac', '%d/%d' % (sum(has), len(has)))] += 1
                # connectivity of the alpha-edge adjacency among tiles with an alpha-edge
                T = [t for t in tiles if alpha_edges(t, al)]
                if T:
                    seen = {0}; todo = [0]
                    while todo:
                        i = todo.pop()
                        for j in range(len(T)):
                            if j not in seen and alpha_edges(T[i], al) & alpha_edges(T[j], al): seen.add(j); todo.append(j)
                    st[('zone', name, 'alpha-edge tiles connected via shared alpha-edges', len(seen) == len(T))] += 1
print('SIZE', SIZE)
for kk in sorted(st, key=str): print(st[kk], kk)
print('JOINT (E0,E1,L0,L1) translation classes rel. shared square:', len(classes))
for kk in sorted(classes, key=lambda q: (q[0], q[1], -classes[q])):
    print('CLASS', kk[0], 'T_L', kk[1], 'n', classes[kk], 'types E0 E1 L0 L1', ex[kk])
pickle.dump({kk: classes[kk] for kk in classes}, open('/home/clio/projects/proofs/code-1126c2-Q597/joint%d.pkl' % SIZE, 'wb'))
