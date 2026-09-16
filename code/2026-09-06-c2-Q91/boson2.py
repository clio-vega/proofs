import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine
from bead import trim, t
from boson import Rstar
print("[R_e^*,R_e] -- is it the scalar [e]_{t^2} = 1+t^2+...+t^{2(e-1)} ?")
for e in (1, 2, 3, 4, 5, 6):
    pred = sum(t**(2*h) for h in range(e))
    ok = 0; tot = 0; bad = []
    for n in range(0, 8):
        for lam in (engine.parts_of(n) if n else ((),)):
            lam = trim(lam)
            a = Rstar(engine.op_R({lam: 1}, e), e)
            b = engine.op_R(Rstar({lam: 1}, e), e)
            c = engine.sub(a, b)
            tot += 1
            want = {lam: sp.expand(pred)} if sp.expand(pred) != 0 else {}
            got = {k: sp.expand(v) for k, v in c.items()}
            if got == want: ok += 1
            else: bad.append((lam, got))
    print(f"  e={e}: {ok}/{tot} match  [e]_(t^2) = {sp.expand(pred)}")
    for bb in bad[:2]: print("     MISMATCH", bb)
