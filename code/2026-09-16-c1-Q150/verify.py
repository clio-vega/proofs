"""
Q150 verification suite.  Checks V1-V5 of proofs/2026-09-16-c1-Q150-nonvanishing.tex.
Run: python3 verify.py
"""
import math, sys
import reduce, direct
from reduce import W, tensor, twobead_pairs, twobead_pairs_strong, compat_pairs, sigmabar

def Csup(d, delta, b):
    """C_{delta,b} = { u in {0,1}^{d-1} : u_delta = b, u_{d-delta} = 1-b }."""
    return frozenset(u for u in W(d-1) if u[delta-1] == b and u[d-delta-1] == 1-b)

def wt(sig):   # sigma-coordinates -> the actual weight W
    return {u: ((-1)**sum(u))*sig[u] for u in sig}

# ---------------------------------------------------------------- V1
def V1(pairs, p=3):
    """Structure theorem: every solution of the reduced system over F_p has the
    form sigma = alpha tau^(x)(e/d), rho = beta tau^(x)((f-e)/d) for tau of rank d."""
    print("V1  structure theorem sigma = alpha tau^(x)(e/d), tau of rank d=gcd(e,f)")
    for (e, f) in pairs:
        a, q = e-1, f-e-1
        d = math.gcd(e, f)
        sols = reduce.solve_Fp(a, q, p)
        bad = 0
        for sig, rho in sols:
            hit = False
            for tv in __import__('itertools').product(range(p), repeat=2**(d-1)):
                if not any(tv): continue
                tau = dict(zip(W(d-1), tv))
                S = tensor(tau, d, e//d); R = tensor(tau, d, (f-e)//d)
                # sigma = alpha S, rho = beta R  for some alpha,beta in F_p^x
                for al in range(1, p):
                    if all((sig[u] - al*S[u]) % p == 0 for u in S):
                        for be in range(1, p):
                            if all((rho[v] - be*R[v]) % p == 0 for v in R):
                                hit = True; break
                    if hit: break
                if hit: break
            if not hit: bad += 1
        print(f"     (e,f)=({e},{f}) d={d}: {len(sols):4d} solutions over F_{p}, "
              f"{bad} NOT of the tau-form")
        assert bad == 0

# ---------------------------------------------------------------- V2
def V2(cases):
    """Counterexamples for d>=3: tau = 1_{C_{1,b}} lifts to a commuting pair with zeros."""
    print("V2  counterexamples at d=gcd(e,f)>=3   (tau = indicator of C_{1,b})")
    for (e, f, b) in cases:
        d = math.gcd(e, f); a, q = e-1, f-e-1
        C = Csup(d, 1, b)
        tau = {u: (1 if u in C else 0) for u in W(d-1)}
        sig = tensor(tau, d, e//d); rho = tensor(tau, d, (f-e)//d)
        n2, nC = reduce.check_system(e, f, sig, rho)
        zs = sum(1 for v in sig.values() if v == 0)
        sgb = sigmabar(sig, rho, a, q)
        zb = sum(1 for v in sgb.values() if v == 0)
        v = "(window too large)"
        if e+f <= 15:
            v = "COMMUTES" if not direct.commutator_nonzero(e, f, wt(sig), wt(sgb),
                                                            Lw=10, cap=1) else "FAILS"
        print(f"     (e,f)=({e},{f}) d={d} b={b}: (2B*_e) viol={n2}, (C) viol={nC}; "
              f"W has {zs}/{2**a} zeros, Wbar has {zb}/{2**(f-1)} zeros; {v}")
        assert n2 == 0 and nC == 0 and zs > 0 and zb > 0

# ---------------------------------------------------------------- V3
def V3(pairs, p=5):
    """d <= 2  =>  every solution is nowhere zero (exhaustive over F_p)."""
    print(f"V3  gcd(e,f)<=2  =>  nonvanishing   (exhaustive over F_{p})")
    for (e, f) in pairs:
        a, q = e-1, f-e-1
        assert math.gcd(e, f) <= 2
        sols = reduce.solve_Fp(a, q, p)
        bad = [1 for s, r in sols if 0 in s.values() or 0 in r.values()]
        print(f"     (e,f)=({e},{f}) d={math.gcd(e,f)}: {len(sols):5d} solutions, "
              f"{len(bad)} with a vanishing weight value")
        assert not bad

# ---------------------------------------------------------------- V4
def V4():
    """Negative control: the lift must FAIL when tau's rank does not divide BOTH
    e and f, and when tau violates the rank-d two-bead condition."""
    print("V4  negative controls (these MUST fail)")
    # (a) tau of rank 3 lifted to e=6 but f=10 (3 does not divide 10)
    tau = {u: (1 if u in Csup(3,1,1) else 0) for u in W(2)}
    sig = tensor(tau, 3, 2)                       # rank 6
    rho = {v: (1 if v[0:2] in Csup(3,1,1) else 0) for v in W(3)}   # ad hoc rank 4
    n2, nC = reduce.check_system(6, 10, sig, rho)
    print(f"     (a) d=2 but tau of rank 3 at (6,10): (2B*) viol={n2}, (C) viol={nC}"
          f"   -> {'control fires' if (n2 or nC) else 'CONTROL DEAD'}")
    assert n2 or nC
    # (b) tau = 1_{(0,0)} at d=3 violates (2B*_3) -- must break (2B*_e) after lifting
    tau2 = {u: (1 if u == (0,0) else 0) for u in W(2)}
    s2 = tensor(tau2, 3, 2); r2 = tensor(tau2, 3, 1)
    n2b, nCb = reduce.check_system(6, 9, s2, r2)
    print(f"     (b) tau=1_(0,0) fails (2B*_3), lifted to (6,9): (2B*) viol={n2b}, "
          f"(C) viol={nCb}   -> {'control fires' if n2b else 'CONTROL DEAD'}")
    assert n2b
    # (c) the commutator itself must be nonzero for (b)
    sgb = sigmabar(s2, r2, 5, 2)
    bad = direct.commutator_nonzero(6, 9, wt(s2), wt(sgb), Lw=10, cap=2)
    print(f"     (c) direct commutator for (b): {len(bad)} nonzero elements"
          f"   -> {'control fires' if bad else 'CONTROL DEAD'}")
    assert bad

# ---------------------------------------------------------------- V5
def V5():
    """(2B*_e) for tau^(x)k  <=>  (2B**_d) for tau when k>=2, (2B*_d) when k=1."""
    print("V5  reduction of (2B*_e) to a rank-d condition on tau")
    for d in (2, 3, 4):
        for k in (1, 2, 3):
            e = k*d
            if 2**(e-1) > 2**8: continue
            tb_e = twobead_pairs(e-1)
            tb_d = twobead_pairs(d-1) if k == 1 else twobead_pairs_strong(d)
            agree = 0; tot = 0; one_way = 0
            import itertools
            for tv in itertools.product(range(3), repeat=2**(d-1)):
                if not any(tv): continue
                tau = dict(zip(W(d-1), tv)); tot += 1
                S = tensor(tau, d, k)
                lhs = not any((S[x]*S[y]-S[z]*S[w]) % 3 for (x,y),(z,w) in tb_e)
                rhs = not any((tau[x]*tau[y]-tau[z]*tau[w]) % 3 for (x,y),(z,w) in tb_d)
                agree += (lhs == rhs)
                if rhs and not lhs: one_way += 1     # would break Lemma 5.1(a)
            name = "(2B*_d)" if k == 1 else "(2B**_d)"
            # (2B**_d) => (2B*_e) always (Lemma 5.1(a)); the CONVERSE fails at d>=4.
            # The known exceptions at d=4, k>=2 are the four singletons 1_{100},1_{110},
            # 1_{001},1_{011} and their nonzero scalar multiples: 8 weights over F_3.
            expected_exc = 8 if (d == 4 and k >= 2) else 0
            print(f"     d={d} k={k} (e={e}): (2B*_e) on tau^(x)k  vs  {name} on tau : "
                  f"{agree}/{tot} over F_3  (expected exceptions: {expected_exc})")
            assert tot - agree == expected_exc
            assert not one_way, "(2B**_d) => (2B*_e) must have NO exception"

if __name__ == "__main__":
    V1([(2,3),(2,4),(3,4),(3,5),(3,6),(4,5),(4,6)])
    V2([(3,6,1),(3,6,0),(3,9,1),(6,9,1),(4,8,1),(5,10,1),(6,12,0),(9,12,1),(6,15,1)])
    V3([(2,3),(2,4),(2,5),(3,4),(3,5),(3,7),(4,5),(4,6),(4,7),(5,6),(4,10)])
    V4()
    V5()
    print("\nall checks passed")
