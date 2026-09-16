import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine
from variants import R_variant, VARIANTS
from bead import trim

DATA = [(e, lam) for e in (2, 3, 4) for n in range(0, 8)
        for lam in (engine.parts_of(n) if n else ((),))]

def profile(w, us):
    return tuple(tuple(sorted((k, sp.srepr(sp.expand(v))) for k, v in
                              R_variant(lam, e, w, us).items())) for e, lam in DATA)

shape = tuple(tuple(sorted((k, sp.srepr(sp.expand(v))) for k, v in
                           engine.op_R({trim(lam): 1}, e).items())) for e, lam in DATA)

print(f"data points: {len(DATA)}  (e in 2,3,4;  |lam| <= 7)\n")
profs = {}
for name, w, us in VARIANTS:
    p = profile(w, us)
    n_ok = sum(1 for i in range(len(DATA)) if p[i] == shape[i])
    profs[name] = p
    print(f"{name}  agrees with shape engine on {n_ok}/{len(DATA)}"
          f"   -> {'PASS' if n_ok == len(DATA) else 'FAIL'}")

print("\n--- pairwise distinctness of the variants themselves ---")
names = list(profs)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        if profs[names[i]] == profs[names[j]]:
            print(f"  IDENTICAL: {names[i].split(' ')[0]} == {names[j].split(' ')[0]}"
                  f"   ({names[i][4:].strip()}  vs  {names[j][4:].strip()})")
