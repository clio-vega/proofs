"""Q594 exchange law test on overlapping B/C pairs (and B/B control).  Same sampler.
dC = W'_C ^ W_C, dB = W'_B ^ W_B.  Tests: E1 dC == dB;  E2 split: addC = W'_C - W_C, remC = W_C - W'_C (and for B);
E3 remC subset W_B ? addC disjoint from W_B? remB subset W_C? ...  printed as a Counter of boolean signatures."""
import sys, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1116c2-Q594')
exec(open('/home/clio/projects/proofs/code-1116c2-Q594/locality.py').read().split('st = Counter()')[0])
from collections import Counter
sig = Counter()
for key, b0 in census():
    for b, DV in states(b0, cap=60):
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind, (tx, ty) in (('BC', 'BC'), ('BB', 'BB')):
                    if kind == 'BB' and not x < y: continue
                    WB, b1 = mv(b, x, tx, DV); WC, b2 = mv(b, y, ty, DV)
                    if not (set(WB) & set(WC)): continue
                    WCp, _ = mv(b1, y, ty, DV); WBp, _ = mv(b2, x, tx, DV)
                    B, C, Bp, Cp = map(set, (WB, WC, WBp, WCp))
                    dC, dB = Cp ^ C, Bp ^ B
                    addC, remC, addB, remB = Cp - C, C - Cp, Bp - B, B - Bp
                    s = (kind, 'dC==dB', dC == dB, 'addC==remB', addC == remB, 'remC==addB', remC == addB,
                         'remC<=B', remC <= B, 'addB<=C', addB <= C, 'addC&B', bool(addC & B), '|add|,|rem|', len(addC), len(remC))
                    sig[s] += 1
for k, v in sorted(sig.items(), key=lambda kv: (kv[0][0], -kv[1])): print(v, k)
