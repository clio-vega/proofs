r"""v_S(Atilde^B_{I,I}) at {w_m = 1}, measured per I.  Type B has h=1, so eps=1 is the
only admissible value.  Needed because prop:short of 2026-10-17-c2 proves the
colength count in type C only; the type-B case is used here as a measured input.
DOMAIN: N=1,2,3; X=B; every m; eps=1; gated-generic rational point.
"""
import sys, os
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code-1017c2-coboundary'))
sys.path.insert(0, os.path.dirname(__file__))
import q484_engine as E
from struct_checks import generic_point, ordt
ok = bad = 0
for N in (1, 2, 3):
    for m in range(1, N + 1):
        w = generic_point(N, m, 1, '3/7', 4000 + N)   # gate uses the type-C root product
        for I in E.subsets(N):
            v = ordt(E.Atilde(N, w, I, I, X='B').subs({E.q: sp.Rational(3, 7)}))
            want = 1 if m in I else 0
            if v == want: ok += 1
            else:
                bad += 1
                print(f"  FAIL N={N} m={m} I={I}: v={v} want={want}")
print(f"type-B boundary diagonal valuation v_S(Atilde^B_II) = [m in I]: {ok} ok, {bad} bad")
