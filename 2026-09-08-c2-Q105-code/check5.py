"""CHECK 5: reproduce Q99 Thm B on a fresh engine, then run it at t -> -t."""
import sympy as sp
from vertex import *
t = sp.symbols('t')

def basis(maxdeg):
    for d in range(maxdeg+1):
        for mu in parts(d):
            yield mu, {mu: sp.Integer(1)}

def thmB(u, m, n, f):
    """LHS: H^u_m H^{1/u}_n f - u^{-1} H^{1/u}_n H^u_m f."""
    cu  = lambda k: 1 - u**k
    cui = lambda k: 1 - u**(-k)
    A = Hmode(Hmode(f, n, cui), m, cu)
    B = Hmode(Hmode(f, m, cu), n, cui)
    return sub(A, smul(1/u, B))

def rhs(u, m, n, f):
    N = m+n
    G = h_alphabet(N, lambda k: 1 + u**k)
    return smul((1-u**-1)*u**-m, pmul(G, f))

for label, u in [("u = t", t), ("u = -t", -t)]:
    ok = bad = nonzero = 0
    for mu, f in basis(3):
        for m in range(-2, 4):
            for n in range(-2, 4):
                L = thmB(u, m, n, f); R = rhs(u, m, n, f)
                D = sub(L, R)
                D = {k: sp.simplify(sp.cancel(v)) for k, v in D.items()}
                D = {k: v for k, v in D.items() if v != 0}
                if D: bad += 1; print("MISMATCH", label, mu, m, n, D)
                else: ok += 1
                if L: nonzero += 1
    print(f"{label}:  {ok} agree, {bad} mismatch  ({nonzero} of them with nonzero defect)")

# read off the two alphabets explicitly
print("\nalphabet read-off (mode n of each object):")
for u, name in [(t, "u=t"), (-t, "u=-t")]:
    twist  = [sp.expand(1-u**k) for k in (1,2,3)]
    defect = [sp.expand(1+u**k) for k in (1,2,3)]
    print(f"  {name:6s} twist alphabet 1-u^n = {twist};  defect alphabet 1+u^n = {defect}")
