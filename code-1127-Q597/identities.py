"""Evaluate (A) on the 8 locally derived case-a middles and (Bk) on the 4 local back windows, directly with delta
(lookstep.delta, the readout change of one step), no census.  Closes grading point (c)."""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1126-Q597')
from lookstep import delta
def comb(terms):
    D = {}
    for sg, s in terms:
        for k, v in delta(s).items():
            if not isinstance(v, tuple): raise ValueError('unpaired triangle in a step: %r' % (v,))
            w = tuple(a + sg * b for a, b in zip(D.get(k, (0, 0)), v))
            if w == (0, 0): D.pop(k, None)
            else: D[k] = w
    return D
MIDS = []; BACKS = []
# --- case a: reuse casea.py's construction, capturing the middles
src = open('/home/clio/projects/proofs/code-1127-Q597/casea.py').read()
src = src.replace("            sq = [t for t in HE", "            MIDS.append((TL, cE0, x, y, mids))\n            sq = [t for t in HE").replace("        print(line)", "        pass")
src = src.replace("print('case-a classes'", "_ = ('case-a classes'")
exec(src)
for TL, cE0, x, y, (Em, ELm, Lm, LEm) in MIDS:
    DE = comb([(1, s) for s in ELm] + [(-1, s) for s in Em]); DL = comb([(1, s) for s in LEm] + [(-1, s) for s in Lm])
    print('case a T_L', TL, 'E0', cE0, 'pair', x, y, '| Delta_E == 0:', not DE, '| Delta_L == 0:', not DL)
# --- case b back windows
src = open('/home/clio/projects/proofs/code-1127-Q597/wbclass.py').read()
src = src.replace("        o = min(p for t in win[0][1] for p in t[1])", "        BACKS.append((TL, x, y, win))\n        o = min(p for t in win[0][1] for p in t[1])")
exec(src.split("print('local windows matched")[0])
for TL, x, y, (e0, e1, f0, f1, l2, m1, m2) in BACKS:
    D = comb([(1, f0), (1, f1), (-1, e0), (-1, e1), (-1, m1), (-1, m2), (1, l2)])
    print('case b back T_L', TL, 'pair', x, y, '| (Bk) == 0:', not D)
print('ALL:', len(MIDS), 'case-a middles,', len(BACKS), 'back windows')
# CONTROL: translate one step (E^L's last) by a unit vector; the identity must then FAIL (instrument live)
def shiftv(s, v):
    f = lambda S: frozenset((t[0], frozenset(tuple(a + b for a, b in zip(p, v)) for p in t[1])) for t in S)
    return (s[0], f(s[1]), f(s[2]))
nz = 0; tot = 0
for TL, x, y, (e0, e1, f0, f1, l2, m1, m2) in BACKS:
    for v in ((1, 0, 0, 0), (0, 0, 1, 0)):
        tot += 1; nz += bool(comb([(1, f0), (1, shiftv(f1, v)), (-1, e0), (-1, e1), (-1, m1), (-1, m2), (1, l2)]))
for TL, cE0, x, y, (Em, ELm, Lm, LEm) in MIDS:
    for v in ((1, 0, 0, 0), (0, 0, 1, 0)):
        tot += 1; nz += bool(comb([(1, s) for s in ELm[:-1]] + [(1, shiftv(ELm[-1], v))] + [(-1, s) for s in Em]))
print('CONTROL (one step translated): identity fails in', nz, 'of', tot)
