"""
Independent verification of Rick's Day 158 memo (2026-09-02), from the DEFINITIONS ONLY.
Clio, 2026-09-03 peer review. No use of Rick's scripts.

Objects (his section 1):
    F(T;u1,u2) = sum_k T^k/k! * A_k(u1) A_k(u2),  A_k(x) = (x+1)_k = (x+1)(x+2)...(x+k)
    E1 = u1+u2, E2 = u1*u2, phi(Y) = 1 + E1 Y + E2 Y^2,  Y = T phi(Y)
    q = 1 - E1 T - 2 T E2 Y,   W = Y/(T q)
u-weight: deg u1 = deg u2 = 1.  ell_k[n] := weight-(n+k) homogeneous part of [T^n] log F.
"""
import sympy as sp
from sympy import Rational as R

u1, u2, T = sp.symbols('u1 u2 T')
N = 14                      # truncation order in T

# ---------- F, as a list of polynomial coefficients c[k] = [T^k] F ----------
c = []
for k in range(N + 1):
    A1 = sp.prod([u1 + j for j in range(1, k + 1)])
    A2 = sp.prod([u2 + j for j in range(1, k + 1)])
    c.append(sp.expand(sp.together(A1 * A2 / sp.factorial(k))))

# ---------- G = F'/F  by series division ----------
fp = [ (m + 1) * c[m + 1] for m in range(N) ]        # [T^m] F'
g = []
for m in range(N):
    s = fp[m] - sum(c[j] * g[m - j] for j in range(1, m + 1))
    g.append(sp.expand(s))

# [T^n] log F = g_{n-1} / n
logF = [sp.Integer(0)] + [sp.expand(g[n - 1] / n) for n in range(1, N + 1)]

# ---------- homogeneous decomposition in (u1,u2) ----------
def homog(poly, d):
    """weight-d homogeneous component of poly in u1,u2"""
    p = sp.Poly(sp.expand(poly), u1, u2)
    out = sp.Integer(0)
    for (a, b), co in p.terms():
        if a + b == d:
            out += co * u1**a * u2**b
    return sp.expand(out)

def topweight(poly):
    if sp.expand(poly) == 0: return None
    p = sp.Poly(sp.expand(poly), u1, u2)
    return max(a + b for (a, b), co in p.terms() if co != 0)

print("=" * 78)
print("CHECK 1 — the weight-labeling caveat.")
print('  His claim: "at u3=0, the top weight of [T^n] log F is n+1 and the sub-top is n."')
print("  Note [T^n]F has top weight 2n, so this asserts a collapse 2n -> n+1.")
print("=" * 78)
print(f"{'n':>3} {'topwt([T^n]F)':>14} {'topwt([T^n]logF)':>18} {'n+1':>5} {'match':>6} "
      f"{'wt-n part nonzero (sub-top)':>28}")
ok_caveat = True
for n in range(1, 11):
    tf = topweight(c[n]); tl = topweight(logF[n])
    subtop_nonzero = sp.expand(homog(logF[n], n)) != 0
    m = (tl == n + 1)
    ok_caveat &= m and subtop_nonzero
    print(f"{n:>3} {tf:>14} {tl:>18} {n+1:>5} {str(m):>6} {str(subtop_nonzero):>28}")
print(f"  => caveat's weight statement holds for n<=10: {ok_caveat}")

# ---------- Y, q, W in E1,E2 then substituted to u1,u2 ----------
E1 = u1 + u2; E2 = u1 * u2
Y = [sp.Integer(0)] * (N + 1)                        # Y = T phi(Y)
for n in range(1, N + 1):
    # iterate: Y_n depends only on lower coefficients
    Ysq = [sum(Y[i] * Y[j - i] for i in range(j + 1)) for j in range(N + 1)]
    Y[n] = sp.expand(( [1] + [0]*N )[n-1] if False else 0)
Y = [sp.Integer(0)] * (N + 1)
for _ in range(N + 2):                                # fixed-point iteration
    Ysq = [sp.expand(sum(Y[i] * Y[j - i] for i in range(j + 1))) for j in range(N + 1)]
    newY = [sp.Integer(0)] * (N + 1)
    for n in range(1, N + 1):
        # [T^n] Y = [T^{n-1}] phi(Y) = delta_{n-1,0} + E1*Y_{n-1} + E2*Ysq_{n-1}
        newY[n] = sp.expand((1 if n == 1 else 0) + E1 * Y[n - 1] + E2 * Ysq[n - 1])
    Y = newY
print("\n  Y_n (n=1..4):", [sp.factor(Y[n]) for n in range(1, 5)],
      " weights:", [topweight(Y[n]) for n in range(1, 5)], "(expect n-1)")

def mul(a, b):
    return [sp.expand(sum(a[i] * b[n - i] for i in range(n + 1))) for n in range(N + 1)]
def inv(a):   # a[0] must be 1
    r = [sp.Integer(0)] * (N + 1); r[0] = sp.Integer(1)
    for n in range(1, N + 1):
        r[n] = sp.expand(-sum(a[j] * r[n - j] for j in range(1, n + 1)))
    return r
def logser(a):  # a[0]=1 ; log a  via (log a)' = a'/a
    ap = [(m + 1) * a[m + 1] for m in range(N)] + [sp.Integer(0)]
    d = mul(ap, inv(a))
    return [sp.Integer(0)] + [sp.expand(d[n - 1] / n) for n in range(1, N + 1)]

# q = 1 - E1 T - 2 T E2 Y   => [T^n] q
qs = [sp.Integer(0)] * (N + 1); qs[0] = sp.Integer(1); qs[1] = sp.expand(-E1)
for n in range(1, N + 1):
    qs[n] = sp.expand(qs[n] - 2 * E2 * Y[n - 1]) if n >= 1 else qs[n]
qs[1] = sp.expand(-E1 - 2 * E2 * Y[0])
for n in range(2, N + 1):
    qs[n] = sp.expand(-2 * E2 * Y[n - 1])
# W = Y/(T q) = (Y/T) / q
YoverT = [sp.expand(Y[n + 1]) for n in range(N)] + [sp.Integer(0)]
W = mul(YoverT, inv(qs))
print("  W_0 =", W[0], "(expect 1);  W_1 =", sp.expand(W[1]))
logW = logser(W)

print("\n" + "=" * 78)
print("CHECK 2 — Prop. A:  T^2 F'' + [(E1+3)T - 1] F' + (1+E1+E2) F = 0")
print("=" * 78)
resid_ok = True
for n in range(0, N - 2):
    # [T^n] of each term
    t1 = (n)*(n-1)*c[n] if n >= 2 else 0                       # T^2 F'' -> n(n-1)c_n
    t2 = (E1 + 3) * (n * c[n]) if n >= 1 else 0                # (E1+3)T F' -> n c_n
    t3 = -(n + 1) * c[n + 1]                                   # -F'
    t4 = (1 + E1 + E2) * c[n]
    r = sp.expand(t1 + t2 + t3 + t4)
    if r != 0:
        resid_ok = False; print(f"   n={n}: residual = {sp.factor(r)}")
print(f"  Prop. A residual identically zero for n=0..{N-3}: {resid_ok}")

print("\n" + "=" * 78)
print("CHECK 3 — Thm 1 (top layer):  [T^n] Xi = E2 * Y_n / n,  Xi = ell_1(log F)")
print("=" * 78)
thm1 = True
for n in range(1, 11):
    lhs = homog(logF[n], n + 1)
    rhs = sp.expand(E2 * Y[n] / n)
    okn = sp.expand(lhs - rhs) == 0
    thm1 &= okn
    if n <= 4 or not okn:
        print(f"   n={n}: ell_1 = {sp.factor(lhs)} | E2 Y_n/n = {sp.factor(rhs)} | {okn}")
print(f"  Thm 1 holds n=1..10: {thm1}")

print("\n" + "=" * 78)
print("CHECK 4 — Thm 2 (sub-top layer):  X^(0) = (1/2) log W,  X^(0) = ell_0(log F)")
print("=" * 78)
thm2 = True
for n in range(1, 11):
    lhs = homog(logF[n], n)
    rhs = sp.expand(logW[n] / 2)
    okn = sp.expand(lhs - rhs) == 0
    thm2 &= okn
    if n <= 4 or not okn:
        print(f"   n={n}: ell_0 = {sp.factor(lhs)} | (1/2)[T^n]logW = {sp.factor(rhs)} | {okn}")
print(f"  Thm 2 holds n=1..10: {thm2}")

import pickle
pickle.dump({'logF':[sp.srepr(x) for x in logF],'logW':[sp.srepr(x) for x in logW],
             'Y':[sp.srepr(x) for x in Y],'W':[sp.srepr(x) for x in W]},
            open('/tmp/day158_state.pkl','wb'))
print("\n[state saved]")
