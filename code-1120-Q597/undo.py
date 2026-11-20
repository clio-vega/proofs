"""PROVE 11-20: Claim U on real journeys.  For every step k and every candidate hexagon S at the state before step k,
classify S by the abstract table of localenum.py (translate so that R's base is 0): good (DL True) / bad (DL False).
Claim U: every bad candidate has tile set == the NEW tiles of step k-1 (the undo of the previous step); never at k=0.
Also: the selected candidate is always good.   Lesion: same statistic with the selected step's own hexagon replaced by n_{k-1} check off."""
import sys, pickle, glob
from collections import Counter
sys.argv = [sys.argv[0], '5', '1', '1']
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590'); from common import unit_dir
from stepset import candidates
import os; os.chdir('/home/clio/projects/proofs/code-1120-Q597')
import localenum as LE
DVEC = {'B': unit_dir((0, 0, 0, 1)), 'C': unit_dir((0, -1, 0, 0))}
def canon(S, R):
    """translate S so that R's lexicographically-least vertex... we use: R's vertex set minus base where base chosen so R=[0,u]x[0,u]"""
    # find base: R = {b, b+U0, b+U1, b+U0+U1}; s0 = pi0 segment; choose u as the difference with sign so that base is pi-min
    pts0 = sorted({P[:2] for P in R[1]}); pts1 = sorted({P[2:] for P in R[1]})
    u0 = tuple(b-a for a, b in zip(*pts0)); u1 = tuple(b-a for a, b in zip(*pts1))
    if u0 == u1: base = pts0[0] + pts1[0]; u = u0
    elif u0 == tuple(-c for c in u1): base = pts0[0] + pts1[1]; u = u0   # s1 = [pts1[1], pts1[1]+u0]
    else: return None, None
    sh = lambda P: tuple(int(a-b) for a, b in zip(P, base))
    return u, frozenset((t[0], frozenset(sh(P) for P in t[1])) for t in S)
TAB = {}
for T in 'BC':
    for u in LE.UNITS:
        Rr, Hs = LE.classify(T, u)
        for typ, H in Hs:
            r = LE.analyse(T, Rr, typ, set(H))
            TAB[(T, u, H)] = (typ, r)
st = Counter()
for fn in sorted(glob.glob('dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC', 'BB2'):
        for TT, M0, J in D[fam]:
            Bd = set(M0); d = DVEC[TT]; prev_new = None
            for k, (typ, g, n) in enumerate(J):
                R = next(t for t in g if t[0] == 'rh')
                for (ty, S, c) in candidates(Bd, R, d):
                    u, cs = canon(S, R)
                    ent = TAB.get((TT, u, cs))
                    if ent is None: st[(fam, 'candidate NOT in abstract table')] += 1; continue
                    if not ent[1]['xrule']: st[(fam, 'candidate fails abstract xrule')] += 1
                    good = ent[1].get('DL')
                    sel = (S == g)
                    st[(fam, 'selected' if sel else 'other', 'good' if good else 'BAD', 'across', ent[1]['across'],
                        'undo of k-1' if (prev_new is not None and S == prev_new) else ('k=0' if k == 0 else 'not undo'))] += 1
                Bd -= g; Bd |= n; prev_new = n
for k in sorted(st, key=str): print(st[k], k)
