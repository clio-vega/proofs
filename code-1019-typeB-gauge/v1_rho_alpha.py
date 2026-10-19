"""V1.  rho^alpha at every positive root of B_N, N=2,3, TYPE-B normalisation,
measured on H_alpha by TWO mechanisms, against prod_{m in D(alpha)} psi_B(w_m)^{-1}.

Positive roots of B_N and their collision data (def:D with h=1):
  e_i - e_j   H={w_i=w_j}      s_a M = M sym-diff {i,j}  (|M cap {i,j}|=1)   D = {}
  e_m         H={w_m=1}        s_a M = M \\ {m}   (m in M)                     D = {m}
  e_a + e_b   H={w_a w_b=1}    s_a M = M \\ {a,b} ({a,b} subset M)            D = {a,b}

Mechanism M1  RESTRICT-then-DIVIDE: substitute the divisor into w, evaluate both
              At's on H_alpha, divide.
Mechanism M2  DIVIDE-then-RESTRICT along a uniformiser t: substitute the
              t-deformation of the divisor, cancel the ratio as a rational
              function of t, then set t=0.  A 0/0 on H_alpha is resolved by
              cancellation in t rather than silently.
No radical uniformiser is used anywhere: all three divisors admit a RATIONAL
parametrisation (w_j=w_i ; w_m=1 ; w_b=1/w_a).
"""
import sympy as sp, itertools
from q484_engine import EB, WX, q, wvars

t = sp.symbols('t')
psiB = lambda y: (1-q*y)/(1-q)

def AtB(N, w, M, I):
    s, r = len(M), len(I)
    if s < r: return sp.Integer(0)
    if s == 0: return sp.cancel(1/EB(N, 0, 0, [], w))
    u = [w[m-1] for m in M]
    return sp.cancel(sp.together(WX(N, s, r, u, w, I, 'B')/EB(N, s, r, u, w)))

def roots(N):
    out = []
    for i, j in itertools.combinations(range(1, N+1), 2):
        out.append(("A  e%d-e%d" % (i, j), ("A", i, j)))
    for m in range(1, N+1):
        out.append(("bdy  e%d" % m, ("B", m, None)))
    for a, b in itertools.combinations(range(1, N+1), 2):
        out.append(("pair e%d+e%d" % (a, b), ("P", a, b)))
    return out

def data(N, kind, a, b):
    """returns (strata test, s_alpha, exact subs, t-deformed subs, D-label-set)"""
    w = wvars(N)
    if kind == "A":
        return (lambda M: len({a, b} & set(M)) == 1,
                lambda M: tuple(sorted(set(M) ^ {a, b})),
                {w[b-1]: w[a-1]}, {w[b-1]: w[a-1]*(1+t)}, ())
    if kind == "B":
        m = a
        return (lambda M: m in M,
                lambda M: tuple(x for x in M if x != m),
                {w[m-1]: sp.Integer(1)}, {w[m-1]: 1+t}, (m,))
    return (lambda M: {a, b} <= set(M),
            lambda M: tuple(x for x in M if x not in (a, b)),
            {w[b-1]: 1/w[a-1]}, {w[b-1]: (1+t)/w[a-1]}, (a, b))

rows = []
for N in (2, 3):
    w = wvars(N)
    for rname, (kind, a, b) in roots(N):
        instrat, salpha, sub0, subt, Dl = data(N, kind, a, b)
        pred = sp.cancel(sp.prod([1/psiB(w[m-1]) for m in Dl]).subs(sub0)) if Dl else sp.Integer(1)
        for ssz in range(0, N+1):
            for M in itertools.combinations(range(1, N+1), ssz):
                if not instrat(M): continue
                sM = salpha(M)
                for r in range(0, N+1):
                    for I in itertools.combinations(range(1, N+1), r):
                        A_M = AtB(N, w, M, I); A_s = AtB(N, w, sM, I)
                        lhs0 = sp.cancel(sp.together(A_M.subs(sub0)))
                        rhs0 = sp.cancel(sp.together(A_s.subs(sub0)))
                        # ---- vacuous stratum: |s_alpha M| < r  -> both must vanish
                        if len(sM) < r:
                            rows.append((N, rname, "VAC |sM|<r", f"M={M},I={I}",
                                         sp.simplify(lhs0) == 0 and sp.simplify(rhs0) == 0, False))
                            continue
                        if sp.simplify(rhs0) == 0:
                            rows.append((N, rname, "vacuous rhs=0", f"M={M},I={I}",
                                         sp.simplify(lhs0) == 0, False)); continue
                        m1 = sp.cancel(sp.together(lhs0/rhs0))
                        okM1 = sp.simplify(sp.cancel(m1 - pred)) == 0
                        # ---- mechanism 2
                        ratio = sp.cancel(sp.together(A_M/A_s))
                        rt = sp.cancel(sp.together(ratio.subs(subt)))
                        m2 = sp.cancel(sp.together(rt.subs(t, 0)))
                        okM2 = sp.simplify(sp.cancel(m2 - pred)) == 0
                        agree = sp.simplify(sp.cancel(m1 - m2)) == 0
                        rows.append((N, rname, "M1&M2", f"M={M},I={I}",
                                     okM1 and okM2 and agree, True))

for N in (2, 3):
    print(f"--- N={N}")
    for rname, _ in roots(N):
        for cat in ("M1&M2", "vacuous rhs=0", "VAC |sM|<r"):
            sub = [x for x in rows if x[0] == N and x[1] == rname and x[2] == cat]
            if not sub: continue
            print(f"  {rname:14s} {cat:14s} {sum(1 for x in sub if x[4])}/{len(sub)}")
            for x in sub:
                if not x[4]: print("     FAIL", x)
disc = [x for x in rows if x[5]]
print(f"\nTOTAL discriminating (rhs != 0, both mechanisms, agreeing): "
      f"{sum(1 for x in disc if x[4])}/{len(disc)}")
vac = [x for x in rows if not x[5]]
print(f"TOTAL vacuous/zero strata (0 = rho*0, reported not scored as content): "
      f"{sum(1 for x in vac if x[4])}/{len(vac)}")

# ---- plants on the headline, partitioned by which block they damage ----
print("\n--- plants on V1 ---")
for plabel, pf in [("type-C shape psi=(1-q w^2)/(1-q)", lambda y: (1-q*y**2)/(1-q)),
                   ("multiplier forced to 1",            None),
                   ("psi not psi^{-1}",                  "inv")]:
    fired = tot = inad = 0
    for N in (2, 3):
        w = wvars(N)
        for rname, (kind, a, b) in roots(N):
            instrat, salpha, sub0, subt, Dl = data(N, kind, a, b)
            true = sp.cancel(sp.prod([1/psiB(w[m-1]) for m in Dl]).subs(sub0)) if Dl else sp.Integer(1)
            if pf is None: bad = sp.Integer(1)
            elif pf == "inv": bad = sp.cancel(sp.prod([psiB(w[m-1]) for m in Dl]).subs(sub0)) if Dl else sp.Integer(1)
            else: bad = sp.cancel(sp.prod([1/pf(w[m-1]) for m in Dl]).subs(sub0)) if Dl else sp.Integer(1)
            if sp.simplify(sp.cancel(bad - true)) == 0: inad += 1; continue
            for ssz in range(0, N+1):
                for M in itertools.combinations(range(1, N+1), ssz):
                    if not instrat(M): continue
                    sM = salpha(M)
                    for r in range(0, len(sM)+1):
                        for I in itertools.combinations(range(1, N+1), r):
                            A_s = AtB(N, w, sM, I).subs(sub0)
                            if sp.simplify(sp.cancel(sp.together(A_s))) == 0: continue
                            A_M = AtB(N, w, M, I).subs(sub0)
                            tot += 1
                            if sp.simplify(sp.cancel(sp.together(A_M - bad*A_s))) != 0: fired += 1
    print(f"  {plabel:36s} fired {fired}/{tot}"
          + (f"   [{inad} root-families INADMISSIBLE: plant == truth there]" if inad else ""))
