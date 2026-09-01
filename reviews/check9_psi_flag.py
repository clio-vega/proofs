"""SECONDARY (a): verify Clio's reading of Rick's Day-151 Psi/Psi^+ FLAG
   ON CLIO'S OWN OBJECTS ONLY. No claim about Rick's s*_mu or frak-s_mu is used.

   Clio's Psi   (section 6): Psi(s_mu)(u)  = det[ (u_i)_{k_j} ] / det[ u_i^{n-j} ],  falling factorial
   Clio's Psi^+ (same, rising): Psi^+(s_mu)(u) = det[ (u_i)^{(k_j)} ] / det[ u_i^{n-j} ]
   with k_j = mu_j + n - j.

   CLAIM TO TEST:  Psi^+(s_mu)(u) = (-1)^{|mu|} Psi(s_mu)(-u).
"""
import sympy as sp
from itertools import count

def partitions_upto(n, maxlen):
    out = [()]
    def rec(rem, cap, pre):
        if len(pre) >= maxlen: return
        for p in range(min(rem, cap), 0, -1):
            out.append(tuple(pre+[p])); rec(rem-p, p, pre+[p])
    rec(n, n, [])
    return sorted(set(out))

def falling(y, k): 
    e = sp.Integer(1)
    for t in range(k): e *= (y - t)
    return e
def rising(y, k):
    e = sp.Integer(1)
    for t in range(k): e *= (y + t)
    return e

def Psi(mu, u, n, ff):
    k = [ (mu[j] if j < len(mu) else 0) + n - 1 - j for j in range(n) ]
    num = sp.Matrix(n, n, lambda i, j: ff(u[i], k[j])).det()
    den = sp.Matrix(n, n, lambda i, j: u[i]**(n-1-j)).det()
    return sp.simplify(sp.cancel(sp.expand(num)/sp.expand(den)))

n = 3
u = sp.symbols('u0 u1 u2')
mus = [m for m in partitions_upto(6, n)]
bad = 0
for mu in mus:
    lhs = Psi(mu, u, n, rising)
    rhs = (-1)**sum(mu) * Psi(mu, tuple(-x for x in u), n, falling)
    if sp.simplify(sp.expand(lhs - rhs)) != 0:
        bad += 1; print("  MISMATCH", mu)
print(f"Psi^+(s_mu)(u) == (-1)^|mu| Psi(s_mu)(-u):  {len(mus)-bad}/{len(mus)} partitions |mu|<=6, n=3, symbolic over QQ")

# negative control: drop the sign, and use rising on both sides
bad2 = sum(1 for mu in mus if sp.simplify(sp.expand(Psi(mu,u,n,rising) - Psi(mu, tuple(-x for x in u), n, falling))) != 0)
bad3 = sum(1 for mu in mus if sp.simplify(sp.expand(Psi(mu,u,n,rising) - (-1)**sum(mu)*Psi(mu, tuple(-x for x in u), n, rising))) != 0)
print(f"  neg control [sign dropped]      : {bad2}/{len(mus)} mismatches  {'DETECTED' if bad2 else '*** MISSED ***'}")
print(f"  neg control [rising both sides] : {bad3}/{len(mus)} mismatches  {'DETECTED' if bad3 else '*** MISSED ***'}")
