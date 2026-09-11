"""(a) verify the TWO-bead matrix-element formula; (b) negative controls."""
import sympy as sp, itertools
from engine import partitions, commutator, maya, unmaya

NM = 12
W  = [sp.Integer(1)] + [sp.Symbol('w%d'%i) for i in range(1,NM)]
WB = [sp.Integer(1)] + [sp.Symbol('W%d'%i) for i in range(1,NM)]

def build(M0, L=60):
    """M0 a finite set of 'window' beads together with Z_{<0}; charge-corrected."""
    M = set(range(-L-40, 0)) | set(M0)
    need = len([x for x in M if x >= 0]) - len([x for x in range(-L-40,0) if x not in M])
    deep = sorted(x for x in M if x < -L)[:0]
    # remove `need` beads far below
    cand = [x for x in range(-L-40, -L)]
    for x in cand[:need]:
        M.discard(x)
    return frozenset(M)

# (a) two-bead formula:  <M'|[Re,Rf]|M> = Wb_Q w_{P-k} - w_P Wb_{Q+k}
bad = tested = 0
for (e,f) in [(2,3),(2,4),(3,4),(3,5)]:
    tab = {e: W, f: WB}; w = lambda g,N: tab[g][N]
    for d in range(1, e):                    # c = b+d,  b = 0  -> k = +1
        b, c = 0, d
        A = list(range(b+1, c)); Bm = list(range(c+1, b+e)); C = list(range(b+e+1, c+f))
        for sa in itertools.chain.from_iterable(itertools.combinations(A,r) for r in range(len(A)+1)):
          for sm in itertools.chain.from_iterable(itertools.combinations(Bm,r) for r in range(len(Bm)+1)):
            for sc in itertools.chain.from_iterable(itertools.combinations(C,r) for r in range(len(C)+1)):
                M = build({b, c} | set(sa) | set(sm) | set(sc))
                L = 60
                lam = unmaya(M, L)
                Mp = frozenset((set(M) - {b, c}) | {b+e, c+f})
                mu = unmaya(Mp, L)
                if sum(mu) != sum(lam) + e + f:  continue
                P = len(sa) + 1 + len(sm); Q = len(sm) + len(sc); k = 1
                want = sp.expand(WB[Q]*W[P-k] - W[P]*WB[Q+k])
                got  = commutator(lam, e, f, w).get(mu, sp.Integer(0))
                tested += 1
                if sp.expand(got - want) != 0:
                    bad += 1
                    if bad < 5: print("MISMATCH", (e,f,d,sa,sm,sc), "got",got,"want",want)
print("two-bead formula: tested %d, mismatches %d" % (tested, bad))

# (b) CONTROL 1: e = f  -> commutator identically zero (solver must find NO constraint)
ctrl = set()
for n in range(9):
    for lam in partitions(n):
        for mu,v in commutator(lam, 3, 3, lambda g,N: W[N]).items(): ctrl.add(v)
print("control e=f=3: nonzero commutator entries =", len(ctrl), "(expect 0)")

# (b) CONTROL 2: the equations must genuinely VARY in w2 (no kernel in the new direction)
w = lambda g,N: W[N]
E13 = {sp.expand(v) for n in range(9) for lam in partitions(n) for v in commutator(lam,1,3,w).values()}
print("[R1,R3] equations:", sorted(E13, key=str))
print("  after imposing w1=-1:", sorted({sp.expand(q.subs(W[1],-1)) for q in E13} - {0}, key=str))
print("  -> w2 still constrained:", any(W[2] in q.free_symbols for q in E13))

# (b) CONTROL 3: R_3^w actually distinguishes tau2 from tau1^2 (the deformation is not a kernel)
from engine import apply_R
print("R_3^w s_(1) =", apply_R({(1,): sp.Integer(1)}, 3, w))
