"""T2. The structural locus, measured ON the loci the Q559 sampler excluded.

Theorem A (to be proved): with w distinct and nonzero, B(T)=0 for EVERY r-subset T
   <=>  n - |S| < r,     S = { a : exists b (possibly b=a) with w_a w_b = 1 }.
Equivalently 2*mu + nu > n - r, mu = #reciprocal pairs a<b, nu = #{a : w_a = +-1}.

prop:recip is the special case mu=1, nu=0, n <= r+1.  The CHAIN case mu>=2 is new.
Both Psi implementations (headsum / alternant) are used: they are independent.
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, B

def S_of(w):
    n = len(w); S = set()
    for a in range(n):
        for b in range(a, n):
            if w[a]*w[b] == 1: S.add(a); S.add(b)
    return S

def mu_nu(w):
    n = len(w)
    mu = sum(1 for a in range(n) for b in range(a+1, n) if w[a]*w[b] == 1)
    nu = sum(1 for a in range(n) if w[a]*w[a] == 1)
    return mu, nu

def all_B_zero(w, r):
    return all(B(w, T) == 0 for T in itertools.combinations(range(n_of(w)), r))
def n_of(w): return len(w)

random.seed(565)
PLAIN = [F(2), F(3), F(5), F(7), F(-2), F(-3), F(4), F(-5), F(9), F(-7)]

def build(n, pairs, ones):
    """pairs: list of (a,b) to make reciprocal; ones: dict a->+-1. Rest generic."""
    w = [None]*n
    for a, v in ones.items(): w[a] = F(v)
    pool = iter([F(2), F(3), F(5), F(7), F(11), F(13), F(17), F(19), F(23)])
    for (a, b) in pairs:
        v = next(pool); w[a] = v; w[b] = F(1, 1)/v
    extra = iter([F(4), F(9), F(25), F(6), F(10), F(14), F(15), F(21), F(22)])
    for a in range(n):
        if w[a] is None: w[a] = next(extra)
    return w

rows = []
# sweep: n from 1..6, r from 1..n, and a catalogue of planted configurations
for n in range(1, 7):
    cfgs = [("generic", [], {})]
    if n >= 2: cfgs.append(("1 recip pair", [(0,1)], {}))
    if n >= 4: cfgs.append(("2 recip pairs (CHAIN)", [(0,1),(2,3)], {}))
    if n >= 6: cfgs.append(("3 recip pairs (CHAIN)", [(0,1),(2,3),(4,5)], {}))
    if n >= 1: cfgs.append(("w=+1", [], {0:1}))
    if n >= 2: cfgs.append(("w=+1,w=-1", [], {0:1,1:-1}))
    if n >= 3: cfgs.append(("1 pair + w=+1", [(0,1)], {2:1}))
    if n >= 5: cfgs.append(("2 pairs + w=+1", [(0,1),(2,3)], {4:1}))
    for name, pairs, ones in cfgs:
        w = build(n, pairs, ones)
        if len(set(w)) != n or any(x == 0 for x in w): continue
        mu, nu = mu_nu(w); S = S_of(w)
        for r in range(1, n+1):
            Bz = all(B(w, T) == 0 for T in itertools.combinations(range(n), r))
            pred = (n - len(S) < r)
            # test Psi itself on three thetas with pairwise distinct m_j
            thetas = [tuple(range(r)), tuple(2*j for j in range(r)), tuple(3-j for j in range(r))]
            psis = []
            for th in thetas:
                m = [th[j]+j for j in range(r)]
                if len(set(m)) != r: continue
                try: psis.append((psi_headsum(w, list(th)), psi_alt(w, list(th))))
                except ZeroDivisionError: pass
            agree = all(a == b for a, b in psis)
            psi_all_zero = all(a == 0 for a, b in psis)
            rows.append((n, r, name, mu, nu, len(S), Bz, pred, psi_all_zero, agree, len(psis)))

print(f"{'n':>2} {'r':>2} {'config':22} {'mu':>2} {'nu':>2} {'|S|':>3} {'allB=0':>7} {'pred':>5} {'Psi=0':>6} {'impl':>5}")
A_ok = A_bad = 0; impl_ok = impl_bad = 0; strict = []
for (n,r,name,mu,nu,sS,Bz,pred,pz,agree,k) in rows:
    flag = "" if Bz == pred else "   <<< THEOREM A FAILS"
    if Bz == pred: A_ok += 1
    else: A_bad += 1
    if agree: impl_ok += 1
    else: impl_bad += 1
    # strictness: Psi=0 but not all B=0  => CANCELLATION point
    if pz and not Bz and k > 0: strict.append((n,r,name))
    print(f"{n:>2} {r:>2} {name:22} {mu:>2} {nu:>2} {sS:>3} {str(Bz):>7} {str(pred):>5} {str(pz):>6} {'ok' if agree else 'MISMATCH':>5}{flag}")
print()
print(f"Theorem A prediction (all B(T)=0  <=>  n-|S| < r):  {A_ok}/{A_ok+A_bad} cells")
print(f"headsum vs alternant agreement:                     {impl_ok}/{impl_ok+impl_bad} cells")
print(f"cells with Psi=0 on all tested theta but some B(T)!=0 (cancellation): {strict}")
