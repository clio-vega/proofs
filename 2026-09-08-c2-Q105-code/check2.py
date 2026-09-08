"""CHECK 2: on the hyperbola s=1/t.  Engine computes with u_f = 1/t directly.
Tests: (a) two-bead sector vanishes identically;
       (b) one-bead element == 0  iff  kappa == 0;
       (c) if kappa != 0, ord_{t=-1} == 1 exactly, with leading coeff (-1)^N * kappa."""
from abacus import *

inv = 1 / t
tot = zero = nz = bad = 0
twobead_nonzero = 0
orders = {}
for n in range(0, 7):
    for lam in parts(n):
        for e in range(1, 5):
            for f in range(1, 5):
                out, M0, L = commutator(lam, e, f, u_e=t, u_f=inv)
                for Mp, val in out.items():
                    val = sp.cancel(sp.together(val))
                    if val == 0:
                        continue
                    d = M0 ^ Mp
                    if len(d) == 4:
                        twobead_nonzero += 1
                        continue
                    (a,) = M0 - Mp
                    k, N = kappa(M0, a, e, f)
                    tot += 1
                    # order of vanishing at t=-1
                    num, den = sp.fraction(val)
                    num = sp.Poly(sp.expand(num), t)
                    o = 0
                    p = num
                    while p.eval(-1) == 0 and p.degree() > 0:
                        o += 1; p = p.diff(t)
                    # leading coeff of the (t+1)-expansion
                    ser = sp.series(val, t, -1, 2).removeO()
                    lead = sp.simplify(sp.expand(ser).coeff(t + 1, 1)) if o == 1 else None
                    orders[o] = orders.get(o, 0) + 1
                    if k == 0:
                        bad += 1; print("BAD: kappa=0 but nonzero", lam, e, f, a, val)
                    # leading coefficient test
                    c1 = sp.limit(sp.diff(val, t), t, -1) if o >= 1 else None
                    if o == 1 and sp.simplify(c1 - (-1)**N * k) != 0:
                        bad += 1; print("BADLEAD", lam, e, f, a, val, c1, (-1)**N*k)
                    nz += 1
# now the converse: kappa==0 => element is zero
missing = 0
for n in range(0, 7):
    for lam in parts(n):
        for e in range(1, 5):
            for f in range(1, 5):
                out, M0, L = commutator(lam, e, f, u_e=t, u_f=inv)
                # every legal (e+f)-move target
                for b, w in moves(M0, e + f):
                    Mp = frozenset((M0 - {b}) | {b + e + f})
                    k, N = kappa(M0, b, e, f)
                    v = sp.cancel(sp.together(out.get(Mp, 0)))
                    if (k == 0) != (v == 0):
                        missing += 1
                        print("CONVERSE FAIL", lam, e, f, b, "kappa=", k, "val=", v)
                    zero += (k == 0)
print("nonzero one-bead on ts=1:", nz, " order histogram:", orders)
print("two-bead nonzero on ts=1:", twobead_nonzero)
print("kappa==0 cases (predicted zero):", zero, " converse failures:", missing, " bad:", bad)
