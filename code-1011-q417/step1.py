"""Q417 STEP 1 (pre-registered, PROVE.md sec 3.1).

At M=3 (k=1, ell=1, h=2k=2, w=2ell+2=4, alpha=(1^{2M})=(1^6)):
does each fibre contain a cylindric tableau ALL of whose horizontal strips theta
have J_open(theta) = empty (equivalently theta weakly decreasing, theta = 1^a 0^{h-a})?

Such a tableau has Psi_open_T = 1 (empty product), hence ord_{t=1} Psi_open_T = 0,
which Q410's theorem proved is NECESSARY for M >= 3 (target c_alpha(1) = 1 != 0).

PRINT THE WITNESS, NOT THE COUNT.  A 0 is indistinguishable from a dead pattern
(PROVE.md sec 4, bullet 1).  So this script prints, for every fibre:
  - the number of tableaux,
  - the per-tableau multiset of step-vectors theta,
  - whether every theta is weakly decreasing,
  - and for the first few tableaux, the actual path, thetas, J_open sets and weight.
It also runs a POSITIVE CONTROL: a parameter point where an all-weakly-decreasing
tableau provably EXISTS, so a reading of "none" is known to be the instrument talking.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
from korff import paths, J_open, J_cyclic, mvec, weight_path, t
import sympy as sp

def weakly_decreasing(th):
    return all(th[i] >= th[i+1] for i in range(len(th)-1))

def thetas(P):
    return [tuple(P[a+1][i]-P[a][i] for i in range(len(P[a]))) for a in range(len(P)-1)]

def cminus(la, k, w):
    m = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[m-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0

def survey(k, ell, alpha, show=6, label=""):
    h = 2*k; w = 2*ell+2
    print(f"\n### {label}  k={k} ell={ell} h={h} w={w} alpha={alpha}")
    allP = paths(h, w, alpha)
    print(f"    cylindric tableaux (lattice paths) in total: {len(allP)}")
    # group by endpoint = fibre
    fib = {}
    for P in allP:
        fib.setdefault(P[-1], []).append(P)
    print(f"    fibres (endpoints): {len(fib)}")
    anyw = []
    for la in sorted(fib):
        c = cminus(la, k, w)
        Ps = fib[la]
        wd = [P for P in Ps if all(weakly_decreasing(th) for th in thetas(P))]
        tag = "  <-- HAS all-weakly-decreasing witness" if wd else ""
        print(f"    endpoint {la}  c^-={c:+d}  |fibre|={len(Ps)}  #all-wd={len(wd)}{tag}")
        for P in wd[:2]:
            print(f"         WITNESS path {P}  thetas {thetas(P)}  "
                  f"Psi_open={weight_path(P,w,'Psi_open')}")
        anyw += wd
    print(f"    TOTAL all-weakly-decreasing tableaux over all fibres: {len(anyw)}")
    print(f"    --- first {show} tableaux, verbatim (theta / J_open / J_cyclic / m at each step) ---")
    for P in allP[:show]:
        th = thetas(P)
        print(f"      path {P}  c^-={cminus(P[-1],k,w):+d}")
        for a in range(len(th)):
            print(f"         step {a}: x={P[a]} m={mvec(P[a],w)} theta={th[a]} "
                  f"J_open={J_open(th[a])} J_cyc={J_cyclic(th[a])} "
                  f"wd={weakly_decreasing(th[a])}")
        print(f"         Psi_open(T) = {weight_path(P,w,'Psi_open')}   "
              f"Psi_cyc(T) = {weight_path(P,w,'Psi')}")
    return anyw

if __name__ == "__main__":
    # ---- THE QUESTION: M=3, the first informative row -------------------------
    survey(1, 1, (1,)*6, label="Q417 TARGET M=3 (ell=1)")
    # ---- the ell=0 regression anchor (PROVE.md sec 2: NOT evidence) ----------
    survey(1, 0, (1,)*4, label="REGRESSION ANCHOR M=2 (ell=0)")
    # ---- POSITIVE CONTROL: an alpha where an all-wd tableau MUST exist -------
    # alpha = (h, h, ..., h): every theta = (1,1,...,1), which IS weakly decreasing.
    survey(1, 1, (2,)*2, show=2, label="POSITIVE CONTROL alpha=(h^2): all theta=(1,1) are wd")
