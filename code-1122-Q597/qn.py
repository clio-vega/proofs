"""PROVE 11-22 Step 1a: the window bound for the local contact enumeration.
N(T, cA) = R_p (canonical class-cA rhombus of T at base 0) plus the 9 non-rh tiles of its three good T-configurations.
Q^N(T, cA) = { om_{T'}(R_p, R) : R a class-cJ rhombus of T' (any translate), some good T'-configuration at R has a tile
               (rhombus included) whose interior meets the interior of a tile of N }.
Reports |Q^N|, whether any witness translate lies on the boundary of the search box (must be NO, else enlarge),
the l1-diameter D, hence the window length bound 2D+1 (Lemma window: om is monotone, non-stationary steps move it
by >= 1 in l1, stationary steps never consecutive).
usage: python3 qn.py BOX"""
import sys, os, itertools
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
BOX = int(sys.argv[1])
for T in 'BC':
    Tp = OTHER[T]
    for cA in (0, 1):
        RA = canon(T, cA); N = [RA] + [t for typ in TYPES for t in step(T, RA, typ)[0] if t[0] != 'rh']
        Q = set(); onbd = 0; nv = 0
        for cJ in (0, 1):
            RJc = canon(Tp, cJ)
            for v in itertools.product(range(-BOX, BOX+1), repeat=4):
                RJ = tr(RJc, v)
                hit = False
                for typ in TYPES:
                    rem = step(Tp, RJ, typ)[0]
                    if any(t == s or overlap(t, s) for t in rem for s in N): hit = True; break
                if not hit: continue
                nv += 1
                Q.add(offs(RA, RJ, Tp))
                if max(abs(x) for x in v) == BOX: onbd += 1
        D = max(sum(abs(a-b) for a, b in zip(q1, q2)) for q1 in Q for q2 in Q)
        lo = tuple(min(q[i] for q in Q) for i in range(4)); hi = tuple(max(q[i] for q in Q) for i in range(4))
        print(T, 'cA', cA, '|N|', len(N), 'witness translates', nv, '|Q^N|', len(Q), 'on box boundary', onbd,
              'l1-diam D', D, 'window <= 2D+1 =', 2*D+1, 'box lo', lo, 'hi', hi, flush=True)
