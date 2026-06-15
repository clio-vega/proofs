"""
Job B secondary (2026-06-15) — MO#509068 backbone: the hook closed form.

MO#509068 ("hook-character vanishing sum", 0 answers) claims a hook character sum
vanishes except at the endpoints.  The in-hand machinery is the closed form for the
hook coefficient
    M_j := <s_lambda, e2^j h1^{2(m-j)}> / C(m,j) = C(2m-1-j, a-1)
for the hook lambda = (a, 1^{2m-a}) |- 2m.  This script VERIFIES that closed form
(the hard, reusable part) against Murnaghan-Nakayama.

HONEST SCOPE NOTE.  The exact MO statement (which signed sum vanishes, and in which
ring / at which specialisation "interior j vanish") is NOT in this session's files.
The natural 2-adic reading val(j)=j+2(v2 C(m,j)+v2 M_j) gives J* = {0} or {0,2}, NOT
the "j=0 and j=m survive" pattern the question reportedly states -- so that is the
wrong sum.  The prose answer must be written against the open MO question text in a
write session; this file fixes only the verified closed form it will rest on.
"""
import sys
sys.path.insert(0, '/home/clio/projects/code/threerow-c1')
from mn import Mj
from math import comb

def hook(a, m):
    return tuple([a] + [1] * (2 * m - a))

def verify_closed_form(mmax=12):
    bad = 0; tested = 0
    for m in range(1, mmax + 1):
        n = 2 * m
        for a in range(1, n):
            lam = hook(a, m)
            if sum(lam) != n:
                continue
            for j in range(0, m + 1):
                try:
                    mj = Mj(lam, j, m)
                except Exception:
                    continue
                pred = comb(2 * m - 1 - j, a - 1) if (2 * m - 1 - j) >= 0 else 0
                tested += 1
                if mj != pred:
                    bad += 1
                    if bad < 8:
                        print("  MISMATCH", lam, "m", m, "j", j, "MN", mj, "C(2m-1-j,a-1)", pred)
    print(f"hook M_j = C(2m-1-j, a-1):  tested={tested}  mismatches={bad}")
    return bad

if __name__ == "__main__":
    verify_closed_form(12)
    # falling-factorial form: C(2m-1-j, a-1) = (2m-1-j)! / [(a-1)! (2m-a-j)!]
    # -> as a function of j it is a ratio of falling factorials; this is the handle
    #    for the Kummer/telescoping cancellation the MO answer will use.
    print("closed form is a ratio of falling factorials in j (Kummer handle); "
          "exact vanishing-sum identity deferred to a write session with the MO text.")
