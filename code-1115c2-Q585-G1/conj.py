"""Brief step 2, strong form: is (iii) on M_1 CONJUGATE to (iii) on M (step k journeys identical as hexagon sequences), or do only the
landings agree?  For each step k <= m - alpha_top compare the hexagon sequence of the k-th (iii) journey on M and on M_1:
identical / share a nonempty common suffix / disjoint ending.  Split by whether phi moved the source cell."""
from rec import *
from readouts import unique_flock
from exchange import run_ii_prefix
from collections import Counter
def run_iii_logged(b):
    u, w = b.nest[0]['u'], b.nest[0]['w']; E_ = tuple(-a for a in u); N_ = tuple(-a for a in w)
    fl = unique_flock(b, 0, E_, N_); key = {t: (v, box_of(b, 0, b.tiles[t][1], E_, N_)[0]) for t, v in fl.items()}
    inv = {t: ij for ij, t in b.nest[0]['pos'].items()}; out = []
    for t in sorted(fl, key=lambda t: key[t]):
        _, ij, H = mig_rec(b, 0, inv[t], u); out.append((inv[t], ij, [s[1] for s in H]))
    return out
st = Counter()
for key, b0 in census():
    M = run_iii_logged(copy.deepcopy(b0)); b1 = copy.deepcopy(b0); run_ii_prefix(b1); M1 = run_iii_logged(b1)
    for (s, l, H), (s1, l1, H1) in zip(M, M1):
        shifted = s != s1; k = 0
        while k < min(len(H), len(H1)) and H[-1-k] == H1[-1-k]: k += 1
        st[('shifted' if shifted else 'same cell', 'identical' if H == H1 else ('common suffix' if k else 'no common suffix'))] += 1
for x in sorted(st): print(x, st[x])
