import q451 as Q, sympy as sp, itertools, sys

def report(h, w, lam, n, kmax=2, modulus=None):
    mod = h+w if modulus is None else modulus
    eng = Q.SymEngine(n)
    terms = {}
    for ks in itertools.product(range(-kmax,kmax+1), repeat=h):
        if sum(ks)!=0: continue
        M = sp.Matrix(h,h, lambda a,b: eng.E(lam[a]-(a+1)+(b+1)+mod*ks[a]))
        d = sp.expand(M.det())
        if d!=0: terms[ks]=d
    tot = sp.expand(sum(terms.values()))
    direct = Q.direct_cyl_schur(eng, h, w, lam, n)
    d_deg = sum(lam)
    print(f"  h={h} w={w} lam={lam} n={n} mod={mod} kmax={kmax}: "
          f"{len(terms)} nonzero k-terms, keys={sorted(terms.keys())}")
    print(f"    CJT sum == direct path enumeration ? {sp.expand(tot-direct)==0}"
          f"   (direct has {len(sp.Poly(direct,*eng.x).monoms()) if direct!=0 else 0} monomials,"
          f" |CRST|={sp.Poly(direct,*eng.x).eval({v:1 for v in eng.x}) if direct!=0 else 0} tableaux)")
    for ks in sorted(terms):
        ex = Q.schur_expand(eng, terms[ks], d_deg)
        pos = all(c > 0 for c in ex.values())
        neg = all(c < 0 for c in ex.values())
        tag = "SCHUR-POSITIVE" if pos else ("SCHUR-NEGATIVE" if neg else "MIXED SIGN")
        print(f"    k={ks}: {tag}  expansion={ {m:int(c) for m,c in ex.items()} }")
    return terms, tot, direct

print("=== h=2, w=2  (mod 4) ===")
for lam in [(1,0),(1,1),(2,0),(2,1),(2,2),(3,1),(3,2)]:
    if lam[0]-lam[-1] > 2: continue
    report(2,2,lam, max(sum(lam),2))
    print()
