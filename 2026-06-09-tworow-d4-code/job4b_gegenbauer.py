"""JOB 4 (symbolic Prop 1) — verify  g_b(m) = C_b^{(-m)}(-(1+i)/2)  in QQ(i)[m].

g_b(m) := [u^b] (1 + s u + u^2)^m,  s = 1+i.
Gegenbauer GF:  (1 - 2 x t + t^2)^{-alpha} = sum_n C_n^{alpha}(x) t^n.
With x = -s/2 = -(1+i)/2, alpha = -m, t = u:
    (1 - 2 x u + u^2)^{-alpha} = (1 + s u + u^2)^m,
so  g_b(m) = C_b^{(-m)}(-(1+i)/2)  EXACTLY.

We verify this as an identity of polynomials in the symbolic parameter m
(not merely numerically), so PROVE can cite it cleanly.
"""
import sympy as sp
from sympy import gegenbauer, symbols, expand, simplify, I, Rational

m, u = symbols('m u')
s = 1 + I

def g_b_series(b):
    """g_b(m) via closed form (verified engine): sum_l (m)_{b-l}/((b-2l)! l!) s^{b-2l}."""
    def falling(x, k):
        r = sp.Integer(1)
        for i in range(k):
            r *= (x - i)
        return r
    return expand(sum(falling(m, b-l)/(sp.factorial(b-2*l)*sp.factorial(l))*s**(b-2*l)
                      for l in range(0, b//2+1) if b-2*l >= 0))

if __name__ == "__main__":
    x0 = -(1 + I)/2
    print("=== Prop 1: g_b(m) = gegenbauer(b, -m, -(1+i)/2) symbolically in m ===")
    allok = True
    for b in range(0, 15):
        lhs = g_b_series(b)
        rhs = expand(gegenbauer(b, -m, x0))
        d = simplify(expand(lhs - rhs))
        ok = (d == 0)
        allok &= ok
        print(f"  b={b:2d}: match={ok}" + ("" if ok else f"   diff={d}"))
    print(f"\nProp 1 identity holds symbolically for all b<=14: {allok}")

    print("\n=== three-term recurrence b*g_b = s(m-b+1)g_{b-1} + (2m-b+2)g_{b-2} ===")
    rok = True
    for b in range(2, 13):
        lhs = b*g_b_series(b)
        rhs = s*(m-b+1)*g_b_series(b-1) + (2*m-b+2)*g_b_series(b-2)
        ok = simplify(expand(lhs-rhs)) == 0
        rok &= ok
        if not ok:
            print(f"  b={b}: FAIL")
    print(f"recurrence holds symbolically for all 2<=b<=12: {rok}")
