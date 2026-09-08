"""CHECK 4: matrix elements of multiplication by h_N[(1+t)X] in the Schur basis.
  <s_nu | h_N[(1+t)X] . s_lam> = sum_{k=0}^N t^{N-k} c_k,  c_k = # pairs of horizontal
  strips lam < rho < nu with |rho/lam| = k.   (Pieri, twice.  No plethysm code.)
Tests: palindromy c_k = c_{N-k}; order histogram; explicit non-divisible witnesses."""
import sympy as sp
from check3 import parts, ordvan
t = sp.symbols('t')

def horiz_strips(lam, k):
    """all mu >= lam with mu/lam a horizontal strip of size k."""
    L = len(lam) + 1
    lam = list(lam) + [0]*(L - len(lam))
    out = []
    def rec(i, cap, rem, acc):
        if i == L:
            if rem == 0: out.append(tuple(x for x in acc if x > 0))
            return
        lo = lam[i]
        hi = min(cap, lam[i] + rem) if i > 0 else lam[i] + rem
        # horizontal strip: mu_i <= lam_{i-1}  and mu_i >= lam_i
        hi = min(hi, lam[i-1] if i > 0 else lam[i]+rem)
        hi = min(hi, cap)
        for v in range(lo, hi+1):
            rec(i+1, v, rem - (v-lam[i]), acc+[v])
    rec(0, 10**6, k, [])
    return out

def mat(lam, nu, N):
    """coefficient polynomial in t."""
    poly = sp.Integer(0)
    for k in range(N+1):
        c = 0
        for rho in horiz_strips(lam, k):
            for nu2 in horiz_strips(rho, N-k):
                if nu2 == nu: c += 1
        poly += c * t**(N-k)
    return sp.expand(poly)

def allparts(n):
    return list(parts(n))

hist = {}; pal_ok = pal_bad = 0; witnesses = []
for N in range(1, 7):
    for dl in range(0, 5):
        for lam in allparts(dl):
            for nu in allparts(dl + N):
                p = mat(lam, nu, N)
                if p == 0: continue
                cs = [sp.Poly(p, t).coeff_monomial(t**j) for j in range(N+1)]
                if cs == cs[::-1]: pal_ok += 1
                else: pal_bad += 1; print("NOT PALINDROMIC", lam, nu, p)
                o = ordvan(p)
                hist.setdefault(N, {}).setdefault(o, 0)
                hist[N][o] += 1
                if o == 0 and len(witnesses) < 4 and N % 2 == 0:
                    witnesses.append((N, lam, nu, sp.factor(p)))
print("palindromic:", pal_ok, " non-palindromic:", pal_bad)
print("\nN : {order at t=-1 : count}")
for N in sorted(hist): print("  ", N, hist[N])
print("\nwitnesses with order 0 (NOT divisible by 1+t):")
for w in witnesses: print("   N=%d  <s_%s | . | s_%s> = %s" % (w[0], w[2], w[1], w[3]))
print("\nlambda=(), nu=(N):")
for N in range(1, 7):
    print("   N=%d :" % N, sp.factor(mat((), (N,), N)), "   ord =", ordvan(mat((), (N,), N)))
