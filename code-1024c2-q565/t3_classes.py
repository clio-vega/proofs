"""T3. Where does the cancellation live?

Group the r-subsets T by (K,J): K = set of negation pairs touched, J = set of unpaired
elements used.  The 2^|K| subsets in a class are one orbit of the group (Z/2)^K acting by
negating individual pairs.  Compute each CLASS SUM and ask which classes vanish on their own.
"""
import itertools
from fractions import Fraction as F
from eng import psi_pair, admissible, terms

def build(v, u):
    """M in index order: v_1,-v_1,...,v_lam,-v_lam, u_1..u_f ; return M and the index data."""
    M = []
    pairidx = []
    for x in v:
        pairidx.append((len(M), len(M)+1)); M += [x, -x]
    uidx = list(range(len(M), len(M)+len(u)))
    M += list(u)
    return M, pairidx, uidx

def classify(T, pairidx, uidx):
    Ts = set(T)
    K = []; eps = []
    for i, (p, q) in enumerate(pairidx):
        inp = p in Ts; inq = q in Ts
        if inp and inq: return None          # A(T)=0, dead term
        if inp: K.append(i); eps.append(+1)
        elif inq: K.append(i); eps.append(-1)
    J = tuple(sorted(j for j, ix in enumerate(uidx) if ix in Ts))
    return tuple(K), J, tuple(eps)

def report(v, u, r, label):
    M, pairidx, uidx = build(v, u)
    lam, f = len(v), len(u)
    assert admissible(M, require_recip_free=True), "not reciprocal-free"
    theta = list(range(r))
    cls = {}
    dead = 0
    for T, B, A in terms(M, theta):
        c = classify(T, pairidx, uidx)
        if c is None:
            assert A == 0, "classify says dead but A != 0"
            dead += 1; continue
        K, J, eps = c
        cls.setdefault((K, J), []).append((eps, B*A))
    tot = psi_pair(M, theta)
    print(f"--- {label}:  lam={lam} f={f} n={len(M)} r={r}   Psi={tot}   (r>=f+2: {r>=f+2})")
    print(f"    dead terms (A=0): {dead}")
    allz = True
    for (K, J), items in sorted(cls.items()):
        s = sum((x for _, x in items), F(0))
        L = len(K)
        spect = f - len(J)
        if s != 0: allz = False
        print(f"    K={K} J={J}  L={L} spectators={spect}  |orbit|={len(items)}  classsum {'= 0' if s==0 else '!= 0  -> '+str(s)}")
    print(f"    every class vanishes on its own: {allz}")
    print()
    return allz

random_ok = []
# f=0, L=4 : the cell the global-negation sign CANNOT explain
random_ok.append(report([F(2),F(3),F(5),F(7)], [], 4, "A  lam=4 f=0 r=4  (C(4,2)=6 even)"))
# f=0, L=3, r=3
random_ok.append(report([F(2),F(3),F(5)], [], 3, "B  lam=3 f=0 r=3  (C(3,2)=3 odd)"))
# f=0, r=2, lam=3 : classes have L=2
random_ok.append(report([F(2),F(3),F(5)], [], 2, "C  lam=3 f=0 r=2"))
# f=1 : spectators appear
random_ok.append(report([F(2),F(3),F(5)], [F(11,3)], 3, "D  lam=3 f=1 r=3  (spectator class L=3,s=1)"))
random_ok.append(report([F(2),F(3),F(5),F(7)], [F(11,3)], 3, "E  lam=4 f=1 r=3"))
random_ok.append(report([F(2),F(3),F(5)], [F(11,3),F(13,5)], 4, "F  lam=3 f=2 r=4"))
print("all cells: every class vanished individually?", all(random_ok))
