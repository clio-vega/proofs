"""M-convexity and M^natural-convexity testers for finite subsets of Z^m, with
calibration against sets whose status I know in advance."""

def exch_M(W):
    """symmetric exchange axiom (M-convex set). Returns (bool, witness)."""
    S = set(W)
    m = len(next(iter(S)))
    for w in S:
        for v in S:
            for i in range(m):
                if w[i] <= v[i]:
                    continue
                ok = False
                for j in range(m):
                    if v[j] <= w[j]:
                        continue
                    a = list(w); a[i] -= 1; a[j] += 1
                    b = list(v); b[i] += 1; b[j] -= 1
                    if tuple(a) in S and tuple(b) in S:
                        ok = True; break
                if not ok:
                    return False, (w, v, i)
    return True, None

def exch_Mnat(W):
    """M^natural exchange axiom. Returns (bool, witness)."""
    S = set(W)
    m = len(next(iter(S)))
    for w in S:
        for v in S:
            for i in range(m):
                if w[i] <= v[i]:
                    continue
                a0 = list(w); a0[i] -= 1
                b0 = list(v); b0[i] += 1
                if tuple(a0) in S and tuple(b0) in S:
                    continue
                ok = False
                for j in range(m):
                    if v[j] <= w[j]:
                        continue
                    a = list(w); a[i] -= 1; a[j] += 1
                    b = list(v); b[i] += 1; b[j] -= 1
                    if tuple(a) in S and tuple(b) in S:
                        ok = True; break
                if not ok:
                    return False, (w, v, i)
    return True, None

def sums(W):
    return sorted({sum(w) for w in W})

if __name__ == "__main__":
    from itertools import product
    print("=== CALIBRATION: sets whose status I know before running ===")
    cases = [
        # (name, set, expected M, expected M-nat)
        ("U(2,3) matroid bases {110,101,011}", [(1,1,0),(1,0,1),(0,1,1)], True, True),
        ("{200,020} (no midpoint)",            [(2,0,0),(0,2,0)],          False, False),
        ("box [0,2]x[0,1]",                    list(product(range(3),range(2))), False, True),
        ("box cap sum=2 in [0,2]^3",           [w for w in product(range(3),repeat=3) if sum(w)==2], True, True),
        ("{(0,0),(1,1)} diagonal pair",        [(0,0),(1,1)],              False, False),
        ("simplex sum<=2, w>=0, m=2",          [w for w in product(range(3),repeat=2) if sum(w)<=2], False, True),
        ("singleton",                          [(3,1,4)],                  True, True),
    ]
    allok = True
    for name, W, eM, eN in cases:
        gM = exch_M(W)[0]; gN = exch_Mnat(W)[0]
        tag = "OK " if (gM==eM and gN==eN) else "MISCALIBRATED"
        if tag != "OK ": allok = False
        print(f"  [{tag}] {name:40s} M={gM} (exp {eM})   M-nat={gN} (exp {eN})")
    print("  calibration", "PASSED" if allok else "FAILED -- do not trust results below")
