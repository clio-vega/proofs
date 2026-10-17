"""AUDIT: engine normalises type B by E^C.  Motegi defines widetilde{WB} = WB/E^B.
Measure (a) the exact ratio E^B/E^C under u = w_M, and (b) the TRUE type-B
collapse multiplier on H_ab, against the paper's claimed rho (= type C's).
Pure sympy; no reliance on the branching formula."""
import sympy as sp
import q484_engine as E
from q484_engine import q, wvars, F, EC, WX

def EB(N, n, r, u, w):
    """Normalization factor E^B_{N,n,r}(u,w), eq (typeBnormalizationfactor) L1361-1367."""
    num = sp.Integer(1)
    for j in range(n):
        for k in range(n):
            num *= (1 - q * u[k] / u[j])          # prod_{j,k} (1 - q u_j^{-1} u_k)
    for j in range(n):
        for k in range(N):
            num *= (1 - q * u[j] * w[k])
    for j in range(n):
        num *= (1 - q * u[j])                      # EXTRA vs E^C
    den = sp.Integer(1)
    for j in range(n):
        for k in range(j + 1, n):
            den *= (1 - q * u[j] * u[k])
    for j in range(n):
        den *= (1 - q * u[j]**2)                   # EXTRA vs E^C  (u-diagonal)
    for j in range(N):
        for k in range(j + 1, N):                  # STRICT j<k  (vs j<=k in E^C)
            den *= (1 - q * w[j] * w[k])
    for k in range(N):
        den *= (1 - q * w[k])                      # EXTRA vs E^C
    return num / den / (1 - q)**(n - r)

def Atilde_B_true(N, w, M, I):
    s, r = len(M), len(I)
    if s < r: return sp.Integer(0)
    n = s
    if n == 0:
        return sp.cancel(sp.Integer(1) / EB(N, 0, 0, [], w))
    u = [w[m-1] for m in M]
    return sp.cancel(sp.together(WX(N, n, r, u, w, I, 'B') / EB(N, n, r, u, w)))

# ---------- (a) the exact ratio E^B/E^C under the specialisation u = w_M ----------
print("=== (a) E^B/E^C at u = w_M, predicted prod_{k in Mc} (1-q w_k^2)/(1-q w_k) ===")
for N in (2, 3):
    w = wvars(N)
    for M in E.subsets(N):
        if not M: continue
        n = len(M); u = [w[m-1] for m in M]
        for r in range(0, n+1):
            got = sp.cancel(sp.simplify(EB(N, n, r, u, w) / EC(N, n, r, u, w)))
            pred = sp.Integer(1)
            for k in range(1, N+1):
                if k not in M:
                    pred *= (1 - q*w[k-1]**2)/(1 - q*w[k-1])
            ok = sp.simplify(got - sp.cancel(pred)) == 0
            print(f"  N={N} M={sorted(M)} r={r}: match={ok}" + ("" if ok else f"  got={got} pred={sp.cancel(pred)}"))
        break   # ratio is r-independent except via (1-q)^{-(n-r)} which cancels; one r suffices per M
