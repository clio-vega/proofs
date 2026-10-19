"""V4.  Lemma lem:typeindep -- the UNNORMALISED adjunction ratio is
type-independent:

      W^B(R u D)/W^B(R)  ==  W^C(R u D)/W^C(R).

Mechanism: W^B's prefactor prod_k (1+u_k^{-1})^{-1} restricted to the D-slots is
prod_y y/(y+1), and line 3c at gamma=-1 contributes prod_y (y+1); the product is
prod_y y, which is line 3c at gamma=0.  So this is a claim about ONE block, and
the plant damages exactly that block.

Two readings: (S) w symbolic, (Nw) w numeric with q symbolic.
"""
import sympy as sp
import itertools
import random
from q484_engine import WX, q


def z(e):
    return sp.simplify(sp.cancel(sp.together(e))) == 0


def wnum(N, seed):
    rng = random.Random(seed)
    out = []
    while len(out) < N:
        v = sp.Rational(rng.randint(2, 29), rng.choice([1, 1, 3, 5, 7]))
        if v in out or v * v == 1 or any(v * o == 1 for o in out):
            continue
        out.append(v)
    return out


xx, yy = sp.Rational(3), sp.Rational(5)
CASES = [("d1 eps=+1", [sp.Integer(1)]),
         ("d1 eps=-1", [sp.Integer(-1)]),
         ("d2 free", [xx, 1 / xx]),
         ("d2 {+1,-1}", [sp.Integer(1), sp.Integer(-1)]),
         ("d3 {-1,x,1/x}", [sp.Integer(-1), xx, 1 / xx]),
         ("d4 two free", [xx, 1 / xx, yy, 1 / yy])]
RSETS = [[sp.Rational(7, 2)], [sp.Rational(7, 2), sp.Rational(9, 4)]]


def run(mode, plantf=None):
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
                        lb = WX(N, s_ + d, r, Rl + D, w, I, 'B')
                        rb = WX(N, s_, r, Rl, w, I, 'B')
                        lc = WX(N, s_ + d, r, Rl + D, w, I, 'C')
                        rc = WX(N, s_, r, Rl, w, I, 'C')
                        if z(rb) or z(rc):
                            continue
                        target = sp.cancel(lc / rc)
                        if plantf is not None:
                            target = target * plantf(D)
                        rows.append((name, N, r, z(sp.cancel(lb / rb) - target)))
    return rows


for mode in ('S', 'Nw'):
    rows = run(mode)
    print(f"--- V4 mechanism ({mode}) ---")
    for tag, _ in CASES:
        sub = [a for a in rows if a[0] == tag]
        if not sub:
            print(f"    {tag:16s} (not reached)")
            continue
        print(f"    {tag:16s} {sum(1 for a in sub if a[3])}/{len(sub)}")
        for a in sub:
            if not a[3]:
                print("       FAIL", a)

print("\n--- V4 plant: prefactor NOT cancelled by line 3c (factor prod_y (y+1)) ---")
rows = run('Nw', plantf=lambda D: sp.prod([(v + 1) for v in D]))
for tag, _ in CASES:
    sub = [a for a in rows if a[0] == tag]
    if not sub:
        continue
    print(f"    {tag:16s} fired {sum(1 for a in sub if not a[3])}/{len(sub)}")
