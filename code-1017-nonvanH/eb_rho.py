"""(b) The TRUE type-B collapse multiplier on H_ab = {w_a w_b = 1}, measured.
Compares three candidates:
   rho_C          = (1-q)^2 / ((1-q w_a^2)(1-q w_b^2))         [paper's claim]
   rho_B_pred     = (1-q)^2 / ((1-q w_a)(1-q w_b))             [my E^B correction]
   rho_eng_B      = Atilde_engine^B ratio                       [what the paper measured]
All restricted to H by w_b -> 1/w_a."""
import sympy as sp, itertools
import q484_engine as E
from q484_engine import q, wvars, WX, EC
from eb_audit import EB, Atilde_B_true

def restrictH(f, w, a, b):
    return sp.cancel(sp.together(f.subs(w[b-1], 1/w[a-1])))

def report(N, X):
    w = wvars(N)
    rows = []
    for a, b in itertools.combinations(range(1, N+1), 2):
        for M in E.subsets(N):
            if not (a in M and b in M): continue
            Mp = tuple(sorted(set(M) - {a, b}))
            for I in E.subsets(N):
                if X == 'B':
                    AM  = Atilde_B_true(N, w, M, I)
                    AMp = Atilde_B_true(N, w, Mp, I)
                else:
                    AM  = E.Atilde(N, w, M, I, X='C')
                    AMp = E.Atilde(N, w, Mp, I, X='C')
                AMH, AMpH = restrictH(AM, w, a, b), restrictH(AMp, w, a, b)
                if AMH == 0 and AMpH == 0:  # degenerate: cannot distinguish any multiplier
                    continue
                rows.append((a, b, tuple(sorted(M)), tuple(sorted(I)), AMH, AMpH))
    x = w[0]*0  # placeholder
    cands = {}
    for a, b in itertools.combinations(range(1, N+1), 2):
        wa = w[a-1]
        cands[(a,b)] = {
          'rho_C':      sp.cancel((1-q)**2 / ((1-q*wa**2)*(1-q*wa**-2))),
          'rho_B_pred': sp.cancel((1-q)**2 / ((1-q*wa)*(1-q/wa))),
        }
    tally = {k: 0 for k in ('rho_C','rho_B_pred')}
    n = 0
    for (a,b,M,I,AMH,AMpH) in rows:
        n += 1
        for name, rho in cands[(a,b)].items():
            if sp.simplify(sp.cancel(AMH - rho*AMpH)) == 0:
                tally[name] += 1
    print(f"  type {X}, N={N}: non-degenerate pairs = {n};  " +
          "  ".join(f"{k}: {v}/{n}" for k, v in tally.items()))
    return rows

print("=== (b) collapse multiplier on H, non-degenerate domain only ===")
for N in (2, 3):
    report(N, 'C')
    report(N, 'B')
