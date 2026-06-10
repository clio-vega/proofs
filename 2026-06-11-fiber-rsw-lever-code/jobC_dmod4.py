"""
Job C structural test: is the RSW hit  G_lambda(zeta_d) = f^lambda(zeta_d)  governed by the
trichotomy branch d == 2 (mod 4) (the 'rich' branch)?  Predict: HIT iff d in {2,6,10,...}.

Tests d in {2,3,4,5,6,8} on all lambda |- n, n<=NMAX, comparing G_lambda(zeta_d) (master id)
against the fake degree f^lambda(zeta_d), up to a zeta_d^c shift.  High-precision numeric.
"""
import sys
import mpmath as mp
import sympy as sp
from collections import defaultdict
from fiber_engine import partitions, G_poly, t as TSYM
from jobC_rsw import fake_degree_poly, num_eval_poly_in
mp.mp.dps = 50

def zeta_num(d): return mp.e**(2j*mp.pi/d)

def main(NMAX):
    hit = defaultdict(int); seen = defaultdict(int); both0 = defaultdict(int)
    for n in range(1, NMAX+1):
        for lam in partitions(n):
            lam = tuple(lam)
            Pt = G_poly(lam)
            fp, q = fake_degree_poly(lam)
            for d in (2, 3, 4, 5, 6, 8):
                zc = zeta_num(d)
                L = num_eval_poly_in(TSYM, Pt, zc)
                R = num_eval_poly_in(q, fp, zc)
                seen[d] += 1
                Lz = abs(complex(L)) < 1e-30
                Rz = abs(complex(R)) < 1e-30
                if Lz and Rz:
                    hit[d] += 1; both0[d] += 1
                elif (not Lz) and (not Rz):
                    for c in range(d):
                        if abs(complex(L) - complex(zc**c)*complex(R)) < 1e-25:
                            hit[d] += 1; break
    print(f"=== G_lambda(zeta_d) == zeta_d^c * f^lambda(zeta_d)  (n<={NMAX}) ===")
    print(f"  prediction: HIT iff d == 2 (mod 4)  [rich branch]")
    print(f"{'d':>3} {'d mod4':>7} {'hit/seen':>12} {'(both-zero)':>12}  rich?")
    for d in (2, 3, 4, 5, 6, 8):
        rich = (d % 4 == 2)
        print(f"{d:>3} {d%4:>7} {f'{hit[d]}/{seen[d]}':>12} {f'{both0[d]}':>12}  {'YES' if rich else 'no'}")

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 9)
