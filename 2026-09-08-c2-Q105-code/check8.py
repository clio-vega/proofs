"""CHECK 8: the closed form.   On s = 1/t,
      <M'|[R_e(t),R_f(1/t)]|M>  =  (-1)^N (-t)^P ( 1 - (-t)^kappa ),
   P = N - 2A - m_f,   kappa = 2(A+B-N) + m_e + m_f.
Tested against the abacus engine (which knows none of this)."""
from abacus import *

inv = 1/t
ok = bad = 0; kvals = {}
for n in range(0, 7):
    for lam in parts(n):
        for e in range(1, 5):
            for f in range(1, 5):
                out, M0, L = commutator(lam, e, f, u_e=t, u_f=inv)
                for b, w in moves(M0, e + f):
                    Mp = frozenset((M0 - {b}) | {b + e + f})
                    N, A, B, me, mf = stats(M0, b, e, f)
                    k = 2*(A + B - N) + me + mf
                    P = N - 2*A - mf
                    pred = sp.cancel((-1)**N * (-t)**P * (1 - (-t)**k))
                    got = sp.cancel(sp.together(out.get(Mp, 0)))
                    if sp.simplify(got - pred) == 0:
                        ok += 1; kvals[k] = kvals.get(k, 0) + 1
                    else:
                        bad += 1; print("FAIL", lam, e, f, b, "got", got, "pred", pred)
print("closed form  Xi = (-1)^N (-t)^P (1 - (-t)^kappa):", ok, "agree,", bad, "mismatch")
print("kappa values seen:", dict(sorted(kvals.items())))
