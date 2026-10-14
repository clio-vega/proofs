"""
Q460 -- direct verification of the MECHANISM of the written proof.

This is a DIFFERENT instrument from the rank computation: it uses no linear
algebra and no field arithmetic at all.  It checks the combinatorial step --
that the explicit separating arcs prescribed in the proof really do make the
pole-order difference non-constant.

Construction (as in the proof): a prefix arcs on the EVEN columns 2,4,...,4a
and a suffix arcs on the ODD columns 1,3,...,4a-1, perfectly interleaved.

For P != P' pick the least i with P_i != P'_i, and set
   u = P_i, u' = P'_i, wlog compare endpoints;
   R  = (c-1, c+1)       a short odd arc straddling the first differing even endpoint c
   R' = (1, 4a-1)        the full-span odd arc
The proof claims  d(R) != d(R')  where  d(S) = [u cap S != 0] - [u' cap S != 0].
"""
import itertools, sys

def matchings(cols):
    cols = sorted(cols)
    if not cols:
        yield (); return
    first, rest = cols[0], cols[1:]
    for k in range(len(rest)):
        for m in matchings(rest[:k] + rest[k+1:]):
            yield ((first, rest[k]),) + m

def ordered_matchings(cols):
    return [p for m in matchings(cols) for p in itertools.permutations(m)]

def overlap(u, v): return not (u[1] < v[0] or v[1] < u[0])

def separating_arcs(u, up, a):
    """The proof's prescription: return (R, R')."""
    if u[0] != up[0]:
        c = min(u[0], up[0])
    else:
        c = min(u[1], up[1])
        # straddle the LARGER second endpoint so the shorter arc is missed
        c = max(u[1], up[1])
    R  = (c-1, c+1)
    Rp = (1, 4*a-1)
    return R, Rp

def check(a):
    Cpre = [2*i for i in range(1, 2*a+1)]        # even: 2,4,...,4a
    Codd = [2*i-1 for i in range(1, 2*a+1)]      # odd : 1,3,...,4a-1
    Ps = ordered_matchings(Cpre)
    tot = good = 0
    bad = []
    for P, Pp in itertools.combinations(Ps, 2):
        tot += 1
        i = next(k for k in range(a) if P[k] != Pp[k])
        u, up = P[i], Pp[i]
        R, Rp = separating_arcs(u, up, a)
        # both must be realizable arcs on the odd column set
        okR  = R[0] in Codd and R[1] in Codd and R[0] < R[1]
        okRp = Rp[0] in Codd and Rp[1] in Codd and Rp[0] < Rp[1]
        d1 = (1 if overlap(u, R)  else 0) - (1 if overlap(up, R)  else 0)
        d2 = (1 if overlap(u, Rp) else 0) - (1 if overlap(up, Rp) else 0)
        if okR and okRp and d1 != d2:
            good += 1
        else:
            bad.append((P, Pp, u, up, R, Rp, okR, okRp, d1, d2))
    print("  a=%d : %d/%d prefix pairs separated by the PRESCRIBED arcs "
          "(|prefixes| = %d = (2a)!/2^a)" % (a, good, tot, len(Ps)))
    for b in bad[:6]:
        print("     unseparated: P_i=%s P'_i=%s  R=%s R'=%s okR=%s okR'=%s d=(%d,%d)"
              % (b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9]))
    return good, tot, len(Ps)

print("MECHANISM CHECK -- prescribed separating arcs, no arithmetic used")
for a in (1, 2, 3, 4):
    check(a)

print("\nFALLBACK CHECK -- does SOME pair of odd arcs always separate?")
def check_exists(a):
    Cpre = [2*i for i in range(1, 2*a+1)]
    Codd = [2*i-1 for i in range(1, 2*a+1)]
    oddarcs = [(r, s) for r in Codd for s in Codd if r < s]
    Ps = ordered_matchings(Cpre)
    tot = good = 0
    for P, Pp in itertools.combinations(Ps, 2):
        tot += 1
        sep = False
        for i in range(a):
            if P[i] == Pp[i]: continue
            u, up = P[i], Pp[i]
            vals = set()
            for R in oddarcs:
                vals.add((1 if overlap(u,R) else 0) - (1 if overlap(up,R) else 0))
            if len(vals) > 1: sep = True; break
        good += sep
    print("  a=%d : %d/%d prefix pairs separated by SOME odd arc" % (a, good, tot))
for a in (1, 2, 3, 4):
    check_exists(a)

print("\nNEGATIVE CONTROL -- an all-disjoint column split must FAIL to separate")
def check_block(a):
    """Prefix on the FIRST 2a columns, suffix on the LAST 2a: every pair is
    case A, so no separation is possible and the capacity must collapse."""
    Cpre = list(range(1, 2*a+1))
    Csuf = list(range(2*a+1, 4*a+1))
    sufarcs = [(r, s) for r in Csuf for s in Csuf if r < s]
    Ps = ordered_matchings(Cpre)
    tot = good = 0
    for P, Pp in itertools.combinations(Ps, 2):
        tot += 1
        sep = False
        for i in range(a):
            if P[i] == Pp[i]: continue
            u, up = P[i], Pp[i]
            vals = {(1 if overlap(u,R) else 0) - (1 if overlap(up,R) else 0)
                    for R in sufarcs}
            if len(vals) > 1: sep = True; break
        good += sep
    print("  a=%d : %d/%d separated under the BLOCK split  -> %s"
          % (a, good, tot, "control FIRED (no separation)" if good == 0
             else "control did not fire"))
for a in (2, 3):
    check_block(a)
