"""Is there a GKM-type congruence at the (Q)-divisors?

Colength budgets measured at N=2 (t3) match these pairings:
  (Q1)  {q w_c w_d = 1}  colength #{M : M n {c,d} = {} }  -> pair M <-> M u {c,d}
  (Q1') {q w_c^2  = 1}   colength #{M : c not in M}       -> pair M <-> M u {c}
  (Q2)  {w_c = q w_e}    colength #{I: e in I, c notin I} -> pair M <-> (M\e) u {c}
Measure the ratio Atil_{M',I}/Atil_{M,I} ON the divisor and ask whether it is
independent of I (and of M).  If yes, the multiplier is a 1-cochain and the whole
Theorem-B machinery transfers.  If no, say exactly how it fails.
"""
import sympy as sp, itertools, sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1018c2-gauge')
from q484_engine import q, wvars, Atilde

def subsets(N):
    return [tuple(s) for k in range(N+1) for s in itertools.combinations(range(1,N+1),k)]

def ratios(N, dname, sub, pairs):
    w = wvars(N)
    print(f"\n### N={N}  D = {dname}")
    for (M, Mp) in pairs:
        row = []
        for I in subsets(N):
            a  = sp.cancel(sp.together(Atilde(N, w, Mp, I, X="C").subs(sub)))
            b  = sp.cancel(sp.together(Atilde(N, w,  M, I, X="C").subs(sub)))
            if b == 0 and a == 0: row.append((I, "0/0")); continue
            if b == 0:            row.append((I, "a/0 !!")); continue
            row.append((I, sp.simplify(sp.cancel(a/b))))
        print(f"  M={M} -> M'={Mp}")
        for I, v in row:
            print(f"      I={str(I):10s}  ratio = {v}")

if __name__ == '__main__':
    N = 2; w = wvars(N)
    # (Q1) {q w1 w2 = 1}: w2 = 1/(q w1).  pairing M <-> M u {1,2}
    ratios(N, "{q w1 w2 = 1}", {w[1]: 1/(q*w[0])}, [((), (1,2))])
    # (Q1') {q w1^2 = 1}: w1 = q^{-1/2}.  pairing M <-> M u {1}
    s = sp.sqrt(1/q)
    ratios(N, "{q w1^2 = 1}", {w[0]: s}, [((), (1,)), ((2,), (1,2))])
    # (Q2) {w1 = q w2}.  pairing M <-> (M \ 1) u {2}  and  M <-> (M\2) u {1}
    ratios(N, "{w1 = q w2}", {w[0]: q*w[1]}, [((1,), (2,)), ((2,), (1,))])
