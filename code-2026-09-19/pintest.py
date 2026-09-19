"""B2 PIN TEST for HMMS 1906.09633.

HMMS state their LR corollary for skew shapes nu/kappa with AT MOST ONE BOX PER COLUMN,
where c^nu_{kappa lambda} = K_{lambda, nu-kappa}.  Question: is there ANY instance of
their LR corollary that is not a Kostka statement?

Claim to test: (i) one-box-per-column <=> kappa_i >= nu_{i+1} for all i;
(ii) on that locus c^nu_{kappa,lambda} = K_{lambda, nu-kappa} always;
(iii) the value depends ONLY on (lambda, nu-kappa) -- nu itself is a redundant label.
"""
from lr import lr_coeff
from lorentzian import kostka
from itertools import combinations

def partitions(n, maxpart=None, maxlen=None):
    if maxpart is None: maxpart = n
    out = []
    def rec(rem, mx, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        if maxlen is not None and len(cur) >= maxlen: return
        for p in range(min(rem, mx), 0, -1):
            rec(rem-p, p, cur+[p])
    rec(n, maxpart, [])
    return out

def onebox_per_column(nu, kap):
    """skew nu/kap has at most one box in each column."""
    n = len(nu)
    kap = tuple(kap) + (0,)*(n-len(kap))
    cols = []
    for i in range(n):
        for j in range(kap[i], nu[i]):
            cols.append(j)
    return len(cols) == len(set(cols))

print("=== (i) one-box-per-column  <=>  kappa_i >= nu_{i+1} ===")
bad = 0; tot = 0
for N in range(1, 9):
    for nu in partitions(N):
        n = len(nu)
        for K in range(0, N+1):
            for kap in partitions(K, maxlen=n) + ([()] if K==0 else []):
                kp = tuple(kap)+(0,)*(n-len(kap))
                if len(kap) > n: continue
                if any(kp[i] > nu[i] for i in range(n)): continue
                tot += 1
                lhs = onebox_per_column(nu, kp)
                rhs = all(kp[i] >= (nu[i+1] if i+1 < n else 0) for i in range(n))
                if lhs != rhs: bad += 1
print("   checked %d skew shapes |nu|<=8 ; mismatches = %d" % (tot, bad))

print("\n=== (ii)+(iii) on the one-box-per-column locus: c^nu_{kappa,lambda} = K_{lambda,nu-kappa} ? ===")
viol = 0; checked = 0; seen = {}
for N in range(1, 9):
    for nu in partitions(N):
        n = len(nu)
        for K in range(0, N):
            for kap in ([()] if K==0 else partitions(K, maxlen=n)):
                kp = tuple(kap)+(0,)*(n-len(kap))
                if len(kap) > n: continue
                if any(kp[i] > nu[i] for i in range(n)): continue
                if not onebox_per_column(nu, kp): continue
                mu = tuple(nu[i]-kp[i] for i in range(n))
                for lam in partitions(N-K):
                    c = lr_coeff(nu, kp, lam)
                    k = kostka(lam, mu)
                    checked += 1
                    if c != k:
                        viol += 1
                        if viol <= 5:
                            print("   MISMATCH nu=%s kap=%s lam=%s : c=%d K=%d"%(nu,kp,lam,c,k))
                    key = (lam, tuple(sorted(mu, reverse=True)))
                    seen.setdefault(key, set()).add(c)
print("   checked %d (nu,kappa,lambda) triples on the locus ; violations = %d" % (checked, viol))
multi = {k:v for k,v in seen.items() if len(v) > 1}
print("   (iii) (lambda, sorted(nu-kappa)) classes with >1 distinct LR value:", len(multi))
