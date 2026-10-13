import q451 as Q, sympy as sp, itertools

def report(h, w, lam, n, kmax=1, modulus=None, verbose=True):
    mod = h+w if modulus is None else modulus
    eng = Q.SymEngine(n)
    terms = {}
    for ks in itertools.product(range(-kmax,kmax+1), repeat=h):
        if sum(ks)!=0: continue
        M = sp.Matrix(h,h, lambda a,b: eng.E(lam[a]-(a+1)+(b+1)+mod*ks[a]))
        d = sp.expand(M.det())
        if d!=0: terms[ks]=d
    tot = sp.expand(sum(terms.values()) if terms else sp.Integer(0))
    direct = Q.direct_cyl_schur(eng, h, w, lam, n)
    ok = sp.expand(tot-direct)==0
    ntab = sp.Poly(direct,*eng.x).eval({v:1 for v in eng.x}) if direct!=0 else 0
    print(f"  h={h} w={w} lam={lam} n={n} mod={mod}: {len(terms)} nonzero k-terms | "
          f"reconstruct==direct? {ok} | |CRST(lam;h,w)|={ntab}")
    if not verbose: return ok
    for ks in sorted(terms):
        ex = Q.schur_expand(eng, terms[ks], sum(lam))
        vals=[int(c) for c in ex.values()]
        tag = "POSITIVE" if all(v>0 for v in vals) else ("NEGATIVE" if all(v<0 for v in vals) else "MIXED")
        print(f"     k={ks}: {tag:8s} {len(ex)} Schur terms  { {m:int(c) for m,c in ex.items()} }")
    return ok

print("=== h=3, w=2 (mod 5) ===")
for lam in [(2,1,0),(2,2,1),(3,2,1)]:
    if lam[0]-lam[-1] > 2: continue
    report(3,2,lam, min(sum(lam),4)); print()

print("=== h=2, w=4 (mod 6) ===")
for lam in [(3,0),(4,1),(5,1)]:
    if lam[0]-lam[-1] > 4: continue
    report(2,4,lam, min(sum(lam),4)); print()

print("=== NEGATIVE CONTROL: wrong modulus h+w+1 (must BREAK) ===")
for (h,w,lam) in [(2,2,(3,1)),(2,2,(3,2)),(3,2,(3,2,1)),(2,4,(4,1))]:
    report(h,w,lam, min(sum(lam),4), modulus=h+w+1, verbose=False)
