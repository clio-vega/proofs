"""CHECK 1: abacus engine vs Q96 Thm 4.1(1)+(2), generic symbolic (t,s)."""
from abacus import *

ok = bad = 0; twobead = 0
for n in range(0, 7):
    for lam in parts(n):
        for e in range(1, 5):
            for f in range(1, 5):
                out, M0, L = commutator(lam, e, f)
                for Mp, val in out.items():
                    d = M0 ^ Mp
                    if len(d) == 2:
                        (a,) = M0 - Mp
                        (a2,) = Mp - M0
                        assert a2 - a == e + f
                        pred = q96_one_bead(M0, a, e, f)
                        if sp.expand(val - pred) == 0: ok += 1
                        else: bad += 1; print("MISMATCH", lam, e, f, a, val, pred)
                    else:
                        twobead += 1
print("one-bead agree:", ok, "mismatch:", bad, " (two-bead entries seen:", twobead, ")")
