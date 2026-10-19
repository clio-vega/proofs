"""V2 (cost-aware rewrite).  The general iota-closed adjunction law in type B:
      At_B(R u D) = (prod_{y in D} phi_B(y)) * At_B(R),    phi_B(y)=(1-q)/(1-q y).

TWO MECHANISMS, deliberately different:
  (S)  w symbolic, q symbolic -- exact identity of rational functions.  Only for
       n = |R|+|D| <= 3, because F sums over n! permutations.
  (Nw) w numeric (pairwise-distinct random rationals, none equal to +-1 and no
       two mutually inverse), q SYMBOLIC -- an exact identity in q at a random
       point of the w-torus.  Reaches n = 6.
A green under (S) AND (Nw) is two readings; where only (Nw) reaches, it is
reported as one reading.
"""
import sympy as sp
import itertools
import random
from q484_engine import EB, WX, q
from wbsafe import AtB as _AtB_safe, z as _zc

phiB = lambda v: (1 - q) / (1 - q * v)
psiB = lambda v: (1 - q * v) / (1 - q)


def z(e):
    return _zc(e)


AtB = _AtB_safe   # cancelled form: valid at u = -1 (see wbsafe.py)


def wnum(N, seed):
    rng = random.Random(seed)
    out = []
    while len(out) < N:
        v = sp.Rational(rng.randint(2, 29), rng.choice([1, 1, 3, 5, 7]))
        if v in out or v * v == 1:
            continue
        if any(v * o == 1 for o in out):
            continue
        out.append(v)
    return out


xx, yy = sp.Rational(3), sp.Rational(5)
CASES = [("d0 empty", []),
         ("d1 eps=+1", [sp.Integer(1)]),
         ("d1 eps=-1", [sp.Integer(-1)]),
         ("d2 free {x,1/x}", [xx, 1 / xx]),
         ("d2 {+1,-1}", [sp.Integer(1), sp.Integer(-1)]),
         ("d3 {-1,x,1/x}", [sp.Integer(-1), xx, 1 / xx]),
         ("d4 {x,1/x,y,1/y}", [xx, 1 / xx, yy, 1 / yy])]
RSETS = [[], [sp.Rational(7, 2)], [sp.Rational(7, 2), sp.Rational(9, 4)]]


def run(mode):
    rows = []
    for N in (2, 3):
        w = list(sp.symbols(f'w1:{N+1}')) if mode == 'S' else wnum(N, 77 + N)
        for name, D in CASES:
            for Rl in RSETS:
                s_, d = len(Rl), len(D)
                cap = 3 if mode == 'S' else 6
                if s_ + d > cap:
                    continue
                for r in range(0, s_ + 1):
                    for I in itertools.combinations(range(1, N + 1), r):
                        base = AtB(N, r, Rl, w, I)
                        lhs = AtB(N, r, Rl + D, w, I)
                        rhs = sp.prod([phiB(v) for v in D]) * base
                        rows.append((name, N, r, f"s={s_},I={I}",
                                     z(lhs - rhs), not z(base)))
    print(f"  --- mechanism ({mode}) ---")
    for tag in [c[0] for c in CASES]:
        sub = [a for a in rows if a[0] == tag]
        if not sub:
            print(f"    {tag:18s} (not reached in this mechanism)")
            continue
        disc = [a for a in sub if a[5]]
        print(f"    {tag:18s} green {sum(1 for a in sub if a[4])}/{len(sub)}"
              f"   DISCRIMINATING(base!=0) {sum(1 for a in disc if a[4])}/{len(disc)}")
        for a in sub:
            if not a[4]:
                print("       FAIL", a)
    return rows


print("=== V2 general iota-closed adjunction, type B", flush=True)
rS = run('S')
rN = run('Nw')

print("\n--- plants (admissibility decided by the arithmetic) ---", flush=True)


def plant(badf):
    for name, D in CASES[1:]:
        true = sp.cancel(sp.prod([phiB(v) for v in D]))
        bad = sp.cancel(badf(D))
        if z(true - bad):
            print(f"    {name:18s} INADMISSIBLE (plant multiplier == true multiplier)")
            continue
        fired = tot = 0
        for N in (2, 3):
            w = wnum(N, 77 + N)
            for Rl in RSETS:
                if len(Rl) + len(D) > 6:
                    continue
                for r in range(0, len(Rl) + 1):
                    for I in itertools.combinations(range(1, N + 1), r):
                        base = AtB(N, r, Rl, w, I)
                        if z(base):
                            continue
                        tot += 1
                        if not z(AtB(N, r, Rl + D, w, I) - bad * base):
                            fired += 1
        print(f"    {name:18s} fired {fired}/{tot}")
    print()


print("  P2a: type-C shape (1-q)/(1-q y^2)", flush=True)
plant(lambda D: sp.prod([(1 - q) / (1 - q * v ** 2) for v in D]))
print("  P2b: multiplier forced to 1", flush=True)
plant(lambda D: sp.Integer(1))
print("  P2c: psi instead of psi^{-1}", flush=True)
plant(lambda D: sp.prod([psiB(v) for v in D]))

# ---- the head-exclusion split: r=0 has NO head slot, so head exclusion is
# vacuous there and the law should hold even where head exclusion fails.
print("\n--- split by r (r=0 has no head slot: head exclusion is vacuous) ---")
for tag in [c[0] for c in CASES]:
    for mode, rows in (('S', rS), ('Nw', rN)):
        sub = [a for a in rows if a[0] == tag]
        if not sub: continue
        z0 = [a for a in sub if a[2] == 0]
        zp = [a for a in sub if a[2] > 0]
        print(f"    {tag:18s} ({mode})  r=0: {sum(1 for a in z0 if a[4])}/{len(z0)}"
              f"   r>=1: {sum(1 for a in zp if a[4])}/{len(zp)}")
