"""
Isolate S(n,c) per the TASK definition:
   S(n,c) := a_c / [ (-1)^{n(lam)} f^{lam0} f^{lam1} binom(w,|lam0|) f^{core} ]
where a_c is the leading coefficient of G_lam(q) in powers of (q+1),
c = floor(|2-core|/2), w = #dominoes = (n-|2-core|)/2.

We compute a_c = leading Taylor coeff of G_lam at q=-1.
Check single-valuedness across shapes lam |- n with same c.
Tabulate, conjecture closed form, validate.
"""
import sys
sys.path.insert(0, '/home/clio/projects/scratch/2026-06-04-sT-2quotient')
import sympy as sp
from sympy import binomial, Rational, symbols, factorial, factorint
from syt_lib import (partitions, syt_list, s_stat,
                     two_core_and_quotient, num_syt, n_of_lambda)
q = symbols('q')

def Gpoly(lam):
    n = sum(lam); P = 0
    for U in syt_list(lam):
        P += q**s_stat(U, n)
    return sp.expand(P)

def leading_qplus1(lam):
    """Return (c_detected, a_c) where a_c is leading coeff in (q+1)."""
    G = Gpoly(lam)
    x = symbols('x')
    sh = sp.expand(G.subs(q, x-1))   # substitute q = x-1, so x = q+1
    P = sp.Poly(sh, x)
    # lowest-degree nonzero coefficient
    coeffs = P.all_coeffs()[::-1]  # ascending
    for k, cf in enumerate(coeffs):
        if cf != 0:
            return k, int(cf)
    return None, 0

NMAX = 14
Stab = {}    # (n,c) -> set of (S value, lam) to detect multivaluedness
Sval = {}    # (n,c) -> the single S value (sympy Rational)
detail = {}  # (n,c) -> example lam

print(f"Computing S(n,c) for n=3..{NMAX} ...")
for n in range(3, NMAX+1):
    for lam in partitions(n):
        core, (q0, q1) = two_core_and_quotient(lam)
        c0 = sum(core); c = c0//2; w = (n - c0)//2
        kdet, ac = leading_qplus1(lam)
        assert kdet == c, f"ord mismatch n={n} lam={lam} kdet={kdet} c={c}"
        f0 = num_syt(q0); f1 = num_syt(q1); fcore = num_syt(core)
        denom = ((-1)**n_of_lambda(lam)) * f0 * f1 * int(binomial(w, sum(q0))) * fcore
        assert denom != 0, f"zero denom n={n} lam={lam}"
        Snc = Rational(ac, denom)
        key = (n, c)
        Stab.setdefault(key, set()).add(Snc)
        detail.setdefault(key, []).append((lam, Snc))
    print(f"  n={n} done")

print("\n=== Single-valuedness check ===")
multi = []
for key in sorted(Stab):
    vals = Stab[key]
    if len(vals) > 1:
        multi.append((key, vals))
        print(f"  *** MULTI n={key[0]} c={key[1]}: {vals}")
    else:
        Sval[key] = next(iter(vals))
if not multi:
    print("  ALL (n,c) single-valued. GOOD.")

print("\n=== S(n,c) table ===")
cmax = max(c for (n,c) in Sval)
print("n\\c " + " ".join(f"{c:>14}" for c in range(cmax+1)))
for n in range(3, NMAX+1):
    row = [f"{n:>3}"]
    for c in range(cmax+1):
        if (n,c) in Sval:
            row.append(f"{str(Sval[(n,c)]):>14}")
        else:
            row.append(f"{'.':>14}")
    print(" ".join(row))

# Save Sval for the fitting stage
import pickle
with open('/home/clio/projects/scratch/2026-06-05-Snc/sval.pkl','wb') as fh:
    pickle.dump({k:(int(v) if v==int(v) else v) for k,v in Sval.items()}, fh)
print("\nSaved sval.pkl")
