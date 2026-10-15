"""Q473: the EVALUATION functional and its (reverse) triangularity.

E_kappa := evaluate at x_i = 1/y_{kappa_i}, i = 1..d,  at L = d = |kappa|.

Mechanism: phi_m(1/y_{m'}) = 0  iff  m > m'   (the j=m' factor (1-x y_{m'}) vanishes),
so a (K,sigma) term survives only if every phi index is <= the cap at its variable.

PREDICTED: E_kappa[G_{kappa'}] = 0 unless kappa' <= kappa, where
   kappa' <= kappa  iff  |kappa'| <= |kappa| and kappa'_j <= kappa_j for j <= |kappa'|;
and E_kappa[G_kappa] != 0 with the closed form
   prod_{i<j}(xi_j - q xi_i)/(xi_j - xi_i) * prod h(xi_i) * prod_{i<j}(1-xi_i xi_j)/(1-q xi_i xi_j)
   * prod_i prod_{j<kappa_i} (1-y_j/y_{kappa_i})/(1-q y_j/y_{kappa_i}),   xi_i = 1/y_{kappa_i}.
"""
import engine as E, sympy as sp, random, sys
from sympy import Rational as R
q, a = E.q, E.a

def env(maxy, seed=0):
    rr = random.Random(seed)
    s = {q: R(5,97), a: R(37,101)}
    for j in range(1, maxy+1):
        s[E.ysym(j)] = R(rr.randint(3,29), 103+j)
    return s

def precedes(kp, k):
    """kp <= k : |kp| <= |k| and kp_j <= k_j for j <= |kp|"""
    if len(kp) > len(k):
        return False
    return all(kp[j] <= k[j] for j in range(len(kp)))

def Ekappa(kappa, kp, s):
    """E_kappa[G_{kp}] at L = |kappa|.  G_{kp}=0 identically if |kp| > L (empty subset sum)."""
    d = len(kappa)
    if len(kp) > d:
        return sp.Integer(0)       # vanishing lemma / empty K-sum
    g = E.G(kp, d)
    sub = dict(s)
    for i, m in enumerate(kappa, start=1):
        sub[E.xsym(i)] = 1/s[E.ysym(m)]
    return sp.cancel(sp.together(g.subs(sub)))

def closed_form_diag(kappa, s):
    d = len(kappa)
    xi = [1/s[E.ysym(m)] for m in kappa]
    r = sp.Integer(1)
    for i in range(d):
        r *= E.h_inf(xi[i]).subs(s) if False else s[a]*(1-xi[i]**2)/(s[a]-xi[i])
    for i in range(d):
        for j in range(i+1, d):
            r *= (1-xi[i]*xi[j])/(1-s[q]*xi[i]*xi[j])
            r *= (xi[j]-s[q]*xi[i])/(xi[j]-xi[i])
    for i, m in enumerate(kappa):
        for jj in range(1, m):
            r *= (1-s[E.ysym(jj)]/s[E.ysym(m)])/(1-s[q]*s[E.ysym(jj)]/s[E.ysym(m)])
    return sp.cancel(r)

maxel, maxcard = int(sys.argv[1]), int(sys.argv[2])
s = env(maxel)
kaps = [k for k in E.subsets_upto(maxel, maxcard)]
cells = 0; ok = 0; nz = 0; zc = 0; bad = []
diag_ok = 0; diag_tot = 0; cf_ok = 0
for kap in kaps:
    if len(kap) == 0:
        continue
    for kp in kaps:
        v = Ekappa(kap, kp, s)
        iszero = (sp.simplify(v) == 0)
        pred = precedes(kp, kap)
        cells += 1
        if iszero != pred:
            ok += 1
            (nz, zc) = (nz + (0 if iszero else 1), zc + (1 if iszero else 0))
        else:
            bad.append((kap, kp, iszero, pred))
    dv = Ekappa(kap, kap, s)
    diag_tot += 1
    if sp.simplify(dv) != 0:
        diag_ok += 1
    if sp.simplify(dv - closed_form_diag(kap, s)) == 0:
        cf_ok += 1
    sys.stdout.flush()
print("EVALUATION triangularity, maxel=%d maxcard=%d" % (maxel, maxcard))
print("  %d/%d cells match prediction   (zeros %d, nonzeros %d)" % (ok, cells, zc, nz))
print("  diagonal nonzero: %d/%d" % (diag_ok, diag_tot))
print("  diagonal matches predicted CLOSED FORM: %d/%d" % (cf_ok, diag_tot))
for b in bad[:12]:
    print("  MISMATCH kappa=%s kappa'=%s iszero=%s predicted_le=%s" % b)
