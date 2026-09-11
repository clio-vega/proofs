"""Independent check of the one-bead matrix-element formula, case by case,
against the engine (which knows nothing about routes B/C)."""
import sympy as sp, itertools
from engine import partitions, commutator, maya, unmaya, apply_R

def build_M(e, f, S, L):
    """b=0 ; M = Z_{<0} u {0} u S, charge-corrected by deleting deep beads."""
    g = e + f
    need = 1 + len(S)                 # charge excess
    deep = list(range(-L, -L + need)) # beads to delete, far below the window
    M = set(range(-L - 40, 0)) | {0} | set(S)
    M -= set(deep)
    return frozenset(x for x in M if x >= -L - 40)

def lam_of(M, L):
    return unmaya(M, L)

NMAX = 12
W  = [sp.Symbol('w%d' % i) for i in range(NMAX)]   # w^{(e)}
WB = [sp.Symbol('W%d' % i) for i in range(NMAX)]   # w^{(f)}

def predict(e, f, u, v, z, me, mf):
    """my claimed formula (1-2mf) w_{Sf} Wb_{Pf} - (1-2me) w_{Pe} Wb_{Se}"""
    Pe = u
    Se = v + mf + z
    Pf = u + me + v
    Sf = z
    return sp.expand((1 - 2*mf)*W[Sf]*WB[Pf] - (1 - 2*me)*W[Pe]*WB[Se])

bad = 0; tested = 0
for (e, f) in [(1,2),(1,3),(2,3),(2,4),(3,4),(2,5),(3,5)]:
    g = e + f
    L = 60
    tab = {e: W, f: WB}
    w = lambda gg, N: tab[gg][N]
    A = list(range(1, e)); B = list(range(e+1, f)); C = list(range(f+1, g))
    for me in (0,1):
        for mf in (0,1):
            for sa in itertools.chain.from_iterable(itertools.combinations(A, r) for r in range(len(A)+1)):
                for sb in itertools.chain.from_iterable(itertools.combinations(B, r) for r in range(len(B)+1)):
                    for sc in itertools.chain.from_iterable(itertools.combinations(C, r) for r in range(len(C)+1)):
                        S = set(sa) | set(sb) | set(sc)
                        if me: S.add(e)
                        if mf: S.add(f)
                        M = build_M(e, f, S, L)
                        lam = lam_of(M, L)
                        # target: bead move 0 -> g
                        Mp = (set(M) - {0}) | {g}
                        mup = lam_of(frozenset(Mp), L)
                        got = commutator(lam, e, f, w).get(mup, sp.Integer(0))
                        want = predict(e, f, len(sa), len(sb), len(sc), me, mf)
                        tested += 1
                        if sp.expand(got - want) != 0:
                            bad += 1
                            if bad < 6:
                                print("MISMATCH", (e,f), "S=",sorted(S), "lam=",lam, "mu=",mup, "got",got,"want",want)
print("one-bead formula: tested %d configurations, %d mismatches" % (tested, bad))
