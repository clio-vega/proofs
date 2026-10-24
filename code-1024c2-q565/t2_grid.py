"""T2. Exact grid over the negation stratum: M = {v_i,-v_i}_{i<=lam} u {u_j}_{j<=f},
reciprocal-free, theta=(0,1,...,r-1).  Which (lam,f,r) force Psi=0?

Several independent draws per cell; a cell is reported ZERO only if EVERY draw is 0.
Draws deliberately avoid prod w = +-1 (the accidental locus found in T1) so that
what we measure is the IDENTICAL-VANISHING of the stratum, not a coincidence.
"""
import itertools, random
from fractions import Fraction as F
from eng import psi_pair, psi_head, admissible, terms

random.seed(565_1024)

def draw(lam, f, tries=400):
    for _ in range(tries):
        v = [F(random.randint(2, 40), random.randint(1, 9)) for _ in range(lam)]
        u = [F(random.choice([-1,1])*random.randint(2, 40), random.randint(1, 9)) for _ in range(f)]
        M = []
        for x in v: M += [x, -x]
        M += u
        if not admissible(M, require_recip_free=True): continue
        e = F(1)
        for x in M: e *= x
        if e == 1 or e == -1: continue      # exclude the accidental locus of T1
        # also need negation pairs to be EXACTLY the lam planted ones
        npairs = sum(1 for a in range(len(M)) for b in range(a+1, len(M)) if M[a] == -M[b])
        if npairs != lam: continue
        return M
    return None

print(f"{'lam':>3} {'f':>2} {'n':>2} {'r':>2} | {'r>=f+2':>6} | {'#draws':>6} {'#zero':>5} | nonzero terms (first draw)")
rows = []
for lam in range(0, 5):
    for f in range(0, 4):
        n = 2*lam+f
        for r in range(1, min(n, 6)+1):
            zeros = 0; nd = 0; nzinfo = ""
            for _ in range(6):
                M = draw(lam, f)
                if M is None: continue
                try:
                    a = psi_pair(M, list(range(r)))
                except ZeroDivisionError:
                    continue
                nd += 1
                if nd == 1:
                    nz = sum(1 for _, B, A in terms(M, list(range(r))) if B*A != 0)
                    tot = sum(1 for _ in terms(M, list(range(r))))
                    nzinfo = f"{nz}/{tot}"
                if a == 0: zeros += 1
            if nd == 0: continue
            pred = "yes" if r >= f+2 else "no"
            mark = ""
            if (zeros == nd) != (r >= f+2): mark = "   <<< DISAGREES WITH conj:transversal"
            print(f"{lam:>3} {f:>2} {n:>2} {r:>2} | {pred:>6} | {nd:>6} {zeros:>5} | {nzinfo}{mark}")
            rows.append((lam,f,n,r,zeros==nd,r>=f+2))
print()
dis = [x for x in rows if x[4] != x[5]]
print(f"cells disagreeing with conj:transversal (<=): {len(dis)} of {len(rows)}")
for d in dis: print("   lam,f,n,r =", d[:4], " measured zero:", d[4], " predicted zero:", d[5])
