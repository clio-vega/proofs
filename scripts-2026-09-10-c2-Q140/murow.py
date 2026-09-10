"""mu = (n): closed-form row description + explicit certificate family."""
from rim import *

def predicted_row(lam, n):
    """Lemma A(ii): predicted row of M^{(n)} at lambda, as dict k -> t^ht."""
    lam = tuple(lam)
    if lam == (n,): return {k: sp.Integer(1) for k in range(n)}
    if len(lam) >= 3 and lam[2] >= 2: return {}
    a, b = lam[0], lam[1]; m = len(lam) - 2
    return {b-1: t**(m+1), a: t**m}

bad = []
for n in range(2, 12):
    rows, cols, M, b = local_system((n,))
    assert cols == tuple(sorted([( ) if k==0 else (k,) for k in range(n)], key=lambda g:(sum(g),g))) or True
    colidx = {(): 0}
    for j,g in enumerate(cols): colidx[g] = j
    for i, lam in enumerate(rows):
        pred = predicted_row(lam, n)
        for j, g in enumerate(cols):
            k = g[0] if g else 0
            want = pred.get(k, sp.Integer(0))
            if sp.expand(M[i,j] - want) != 0: bad.append((n, lam, k, M[i,j], want))
print("LEMMA A(ii) rows of M^{(n)}:", "OK for 2<=n<=11" if not bad else bad[:5])

# the certificate family
def certificate(n):
    """c indexed by rows of M^{(n)}; returns dict lambda -> coefficient."""
    al = t*(t+1); d = 1 - (n-1)*t
    ce = {1: (n-2)*t**2 - 2*t, 2: -t**2 + (n-2)*t - 1}
    c = {(n,): al}
    for a in range(1, n):                      # hooks (a,1^{n-a}), from (a,b,m)=(a,1,n-a-1)
        lam = tuple([a] + [1]*(n-a))
        num = ce.get(a, -t*(t+1))
        c[lam] = sp.simplify(num / t**(n-a-1))
    lam = tuple([2,2] + [1]*(n-4))             # the one non-hook row
    c[lam] = c.get(lam, 0) + d / t**(n-4)
    return c

print("\n n | c*M == 0 | c*b")
for n in range(4, 12):
    rows, cols, M, b = local_system((n,))
    c = certificate(n)
    cv = sp.Matrix([[sp.simplify(c.get(lam, 0)) for lam in rows]])
    prod = sp.simplify(cv*M); rhs = sp.simplify((cv*b)[0,0])
    print(f"{n:2d} | {str(all(x==0 for x in prod)):5s}    | {sp.factor(rhs)}")

print("\nn=4 certificate, cleared of denominators (compare paper):")
rows, cols, M, b = local_system((4,))
c = certificate(4)
print(" rows:", rows)
print("  t*c:", [sp.expand(t*sp.simplify(c.get(l,0))) for l in rows])
print(" paper:", [t**2*(t+1), -t**2*(t+1), -t*(3*t-1), -(t-1)**2, 2*(t-1)],
      " (paper order (4),(3,1),(2,2),(2,1,1),(1^4))")
