"""boundbox.py with a reference-pair parameter.  CONTROLS: (i) shpast's common R* (class-0 b rhombus): must be non-empty and
contain the 21 feasible tiles; (ii) lesion (R^A_p, R^J_{c0+3}): non-empty.  Also: which constraint empties the case-a system --
the lambda-part alone (a0 or a1 interval empty for every (s, q, q'))?   usage: python3 boundbox2.py"""
import sys
src = open('boundbox.py').read().split('INF = 10 ** 9')[0]
exec(src)
INF = 10 ** 9
def analyse(T, RA, RJ):
    Tp = OTHER[T]; b = base(RA); nonempty = 0; lamdead = 0; cells = 0; box = [[INF, -INF] for _ in range(4)]
    for s in shapes:
        for q in QS[T]:
            for qp in QS[Tp]:
                cells += 1
                rhs = lambda l, R, qq: qq + mn(l, R) - mn(l, s) - dot(l, b)
                lo = [-INF] * 4; hi = [INF] * 4
                lo[0] = rhs(LAM[0], RA, q[2]); hi[0] = rhs(LAM[0], RJ, qp[2]); lo[2] = rhs(LAM[1], RA, q[3]); hi[2] = rhs(LAM[1], RJ, qp[3])
                if lo[0] > hi[0] or lo[2] > hi[2]: lamdead += 1; continue
                (lA0, lA1), (lJ0, lJ1) = LF[T], LF[Tp]
                cons = [(lA0, '<=', rhs(lA0, RA, q[0])), (lA1, '<=', rhs(lA1, RA, q[1])), (lJ0, '>=', rhs(lJ0, RJ, qp[0])), (lJ1, '>=', rhs(lJ1, RJ, qp[1]))]
                ok = True
                for idx, other in ((1, 0), (3, 2)):
                    l_, h_ = -INF, INF
                    for l, op, c in cons:
                        if l[idx] == 0: continue
                        k = l[other]; sgn = l[idx]
                        if op == '<=':
                            ub = c - (k * lo[other] if k > 0 else k * hi[other])
                            if sgn > 0: h_ = min(h_, ub)
                            else: l_ = max(l_, -ub)
                        else:
                            lb = c - (k * hi[other] if k > 0 else k * lo[other])
                            if sgn > 0: l_ = max(l_, lb)
                            else: h_ = min(h_, -lb)
                    lo[idx], hi[idx] = l_, h_
                    if l_ > h_: ok = False
                if not ok: continue
                nonempty += 1
                for i in range(4): box[i][0] = min(box[i][0], lo[i]); box[i][1] = max(box[i][1], hi[i])
    return cells, lamdead, nonempty, [tuple(int(y) for y in x) for x in box] if nonempty else None
Rstar = CF['B'][(k0, 'b')][0]
for T in 'BC': print('CONTROL shpast R*', T, analyse(T, Rstar, Rstar))
for (T, t0, t1, off) in CASEA:
    Tp = OTHER[T]
    RJ0, RAp, RAp1, RAp2, Hmid = place(T, t0, off)
    H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1]); H2, RJ3, C2 = step(Tp, RJ2, 'b')
    print(T, t0, t1, 'case a (R_p+2, R^J_c0): (cells, killed by lambda alone, non-empty, box)', analyse(T, RAp2, RJ0))
    print('   LESION (R_p, R^J_c0+3):', analyse(T, RAp, RJ3))
    la = lambda R: (min(P[0] for P in R[1]), min(P[2] for P in R[1]))
    print('   lambda(R_p+2) - lambda(R^J_c0) =', tuple(int(x - y) for x, y in zip(la(RAp2), la(RJ0))),
          '| min over Q_T lambda-part', tuple(min(q[i] for q in QS[T]) for i in (2, 3)), 'max over Q_T\' lambda-part', tuple(max(q[i] for q in QS[Tp]) for i in (2, 3)))
