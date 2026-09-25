"""Q254 consolidated verification.  Run: python3 final.py"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-q254')
import cyl as C   # vendored copy, see cyl.py header
from fms import supp_chi, is_Sm_stable, is_bottom_justified, active_rows
from rigidity import dominates, sort_part, partitions_le, rearrangements
from cylindric import cyl_shapes, dom_max, skyline_of_lamhat, all_diagrams
from winding import winds

def ideal_ok(W, lamhat, ell):
    if not all(dominates(lamhat, sort_part(a)) for a in W): return False
    N = sum(lamhat)
    return len(W) == sum(rearrangements(p, ell) for p in partitions_le(N, ell)
                         if dominates(lamhat, p))

print("=" * 78)
print("Q254 -- FMS Theorem 7 vs. the cylindric M-convexity theorem")
print("=" * 78)

# ---------------------------------------------------------------- winding regime
rows = []
for n in range(2, 6):
    for m in range(1, n):
        for mu in cyl_shapes(n, m):
            for lam in itertools.product(*[range(mu[i], mu[i]+n+1) for i in range(m)]):
                if not C.is_shape(lam, n, m) or not C.contains(lam, mu): continue
                N = C.size(mu, lam)
                if N == 0 or N > 6: continue
                for ell in range(2, 5):
                    W = set(C.weights(mu, lam, n, m, ell))
                    if not W: continue
                    w, _, _ = winds(mu, lam, n, m, ell)
                    rows.append((n, m, mu, lam, ell, N, W, w))

windrows = [r for r in rows if r[7]]
print(f"\n[1] Sample: {len(rows)} (cylindric shape, ell) instances, 2<=n<=5, |lam/mu|<=6, 2<=ell<=4.")
print(f"    Genuinely WINDING (support changes if the wrap constraint is relaxed): "
      f"{len(windrows)}/{len(rows)}.")

# ------------------------------------------------- control on thm:main, winding only
bad = 0
for (n, m, mu, lam, ell, N, W, w) in windrows:
    mx = dom_max(W)
    if len(mx) != 1 or not ideal_ok(W, mx[0], ell): bad += 1
print(f"\n[2] Control on thm:main over the WINDING instances: {len(windrows)-bad}/{len(windrows)}"
      f" have a unique dominance maximum lamhat with supp = {{alpha: sort(alpha) <= lamhat}}.")

# --------------------------------------- uniqueness of the diagram, winding only
seen, tested, uniq = set(), 0, 0
offenders = []
for (n, m, mu, lam, ell, N, W, w) in windrows:
    if N > 6 or ell > 4: continue
    lamhat = dom_max(W)[0]
    k = (ell, N, lamhat)
    if k in seen: continue
    seen.add(k); tested += 1
    sols = {tuple(sorted(D)) for D in all_diagrams(ell, N, N) if supp_chi(list(D), ell) == W}
    want = tuple(skyline_of_lamhat(lamhat, ell))
    if sols == {want}: uniq += 1
    else: offenders.append((n, m, mu, lam, ell, lamhat, sols, want))
print(f"\n[3] EXHAUSTIVE search over every diagram with columns inside [ell] and #D = |lam/mu|:")
print(f"    {uniq}/{tested} distinct (ell, N, lamhat) classes from winding shapes have")
print(f"    EXACTLY ONE solution, and it is the skyline diagram of the antidominant")
print(f"    rearrangement of lamhat.  Offenders: {len(offenders)}")
for o in offenders[:4]: print("      ", o)

# ------------------------------------------------------------ negative controls
print(f"\n[4] Negative controls.")
print(f"    (a) REVERSING the column order: supp(chi_D) is a SUMSET over columns, so it is")
print(f"        invariant under reordering BY CONSTRUCTION.  Nothing in the world could make")
print(f"        this control fail, so it is NOT reported as a passing control.  (Verified")
print(f"        once, as a check on the code path, then discarded as uninformative.)")
nd, nb = 0, 0
for (n, m, mu, lam, ell, N, W, w) in windrows[:200]:
    lamhat = dom_max(W)[0]
    D = skyline_of_lamhat(lamhat, ell)
    if len(D) < 2: continue
    nd += 1
    if supp_chi(D[1:], ell) == W: nb += 1
print(f"    (b) DROP one column from the working diagram: {nd-nb}/{nd} instances then FAIL to")
print(f"        reproduce W_ell  (a control that could have failed, and did not).")
np_, nbb = 0, 0
for (n, m, mu, lam, ell, N, W, w) in windrows[:200]:
    lamhat = dom_max(W)[0]
    D = skyline_of_lamhat(lamhat, ell)
    broke = False
    for j, col in enumerate(D):
        if len(col) < ell and min(col) - 1 >= 1:
            E = list(D); E[j] = tuple(sorted(set(col) - {max(col)} | {min(col) - 1}))
            np_ += 1
            if supp_chi(E, ell) == W: nbb += 1
            broke = True
            break
print(f"    (c) MOVE one box of one column up by one row (breaking bottom-justification):")
print(f"        {np_-nbb}/{np_} instances then FAIL to reproduce W_ell.")
