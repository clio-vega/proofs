import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/reviews/2026-09-08-selfreview-code')
from ribbon import border_strips, parts
t, s = sp.symbols('t s')

def conj(lam):
    if not lam: return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))

def basis(N): return sorted(parts(N))

def Rmat(e, N, u):
    """matrix of R_e(u): Lambda_N -> Lambda_{N+e}, in Schur bases."""
    src, tgt = basis(N), basis(N+e)
    M = sp.zeros(len(tgt), len(src))
    for j, lam in enumerate(src):
        for (mu, h) in border_strips(lam, e):
            M[tgt.index(mu), j] += u**h
    return M

def Omega(N):
    """matrix of omega on Lambda_N in the Schur basis."""
    b = basis(N)
    M = sp.zeros(len(b), len(b))
    for j, lam in enumerate(b):
        M[b.index(conj(lam)), j] = 1
    return M

print("Testing   omega R_e(t) omega  ==  t^(e-1) R_e(1/t)   as matrices")
allok = True
for e in (1,2,3,4,5):
    for N in range(0, 9):
        L = Rmat(e, N, t)
        lhs = Omega(N+e) * L * Omega(N)
        rhs = sp.expand(t**(e-1) * Rmat(e, N, 1/t))
        D = sp.simplify(sp.expand(lhs - rhs))
        ok = all(x == 0 for x in D)
        allok &= ok
        if not ok:
            print(f"  FAIL e={e} N={N}"); print(D)
    print(f"  e={e}: OK for N=0..8   (dim {Rmat(e,8,t).shape})")
print("ALL OK" if allok else "FAILURES PRESENT")

# ---------- the two-parameter commutator on the hyperbola ----------
print()
print("Now:  C_{e,f}(t,s) = R_e(t)R_f(s) - R_f(s)R_e(t)  on Lambda_N -> Lambda_{N+e+f}")
print("Predicted:  omega C(t,s) omega = t^(e-1) s^(f-1) C(1/t, 1/s)")
for (e,f) in [(2,2),(2,3),(3,3),(2,4),(3,4)]:
    for N in (0,1,2,3,4):
        A = Rmat(f, N, s); B = Rmat(e, N+f, t)
        A2 = Rmat(e, N, t); B2 = Rmat(f, N+e, s)
        C = B*A - B2*A2
        Ai = Rmat(f, N, 1/s); Bi = Rmat(e, N+f, 1/t)
        Ai2 = Rmat(e, N, 1/t); Bi2 = Rmat(f, N+e, 1/s)
        Ci = Bi*Ai - Bi2*Ai2
        lhs = Omega(N+e+f) * C * Omega(N)
        rhs = sp.expand(t**(e-1)*s**(f-1) * Ci)
        D = sp.simplify(sp.expand(lhs - rhs))
        if not all(x == 0 for x in D):
            print(f"  FAIL (e,f)=({e},{f}) N={N}"); break
    else:
        print(f"  (e,f)=({e},{f}): OK for N=0..4")
