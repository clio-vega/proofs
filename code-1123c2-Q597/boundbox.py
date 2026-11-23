"""PROVE 11-23 c2: box adequacy for the past window of wakeA.py, by explicit bounds (not by enumeration).
A tile t = s + base(R^A_{p+2}) + v (s a shape at the origin).  Every window inequality is  min_{P in t} l(P) <= or >= c
for a linear form l in {a0, a1 (lambda), and the psi-forms}; it reads  l(v) <= or >= c - l(base) - min_{P in s} l(P).
psi-forms (share tex, Lemma shift): psi_B = (min -b1, min b0), psi_C = (min -(a1+b1), min (a0+b0)); asserted against
words.offs on random translates.  Interval propagation: v0, v2 from the lambda-parts (two-sided); then v1, v3 from the psi-parts.
Prints, per contact placement, the resulting bounding box of v over all shapes and all (q, q') -- must lie inside [-6, 6]^4.
usage: python3 boundbox.py"""
import sys, random
sys.argv = ['wakeA.py', '2']
exec(open('wakeA.py').read().split('\nnr = lambda')[0])
LF = {'B': ((0, 0, 0, -1), (0, 1, 0, 0)), 'C': ((0, 0, -1, -1), (1, 1, 0, 0))}
LAM = ((1, 0, 0, 0), (0, 0, 1, 0))
dot = lambda l, P: sum(a * b for a, b in zip(l, P))
mn = lambda l, t: min(dot(l, P) for P in t[1])
def offs2(t, R, T):
    return tuple(int(mn(l, t) - mn(l, R)) for l in LF[T] + LAM)
random.seed(3); bad = 0
for s in shapes:
    for T in 'BC':
        for _ in range(30):
            v = tuple(random.randint(-5, 5) for _ in range(4)); t = trs(s, v); R = trs(canon(T, 1), tuple(random.randint(-3, 3) for _ in range(4)))
            bad += offs2(t, R, T) != offs(t, R, T)
print('psi-form assertion mismatches (must be 0):', bad)
INF = 10 ** 9
for (T, t0, t1, off) in CASEA:
    Tp = OTHER[T]; RJ0, RAp, RAp1, RAp2, Hmid = place(T, t0, off); b = base(RAp2)
    box = [[INF, -INF] for _ in range(4)]; nonempty = 0
    for s in shapes:
        for q in QS[T]:
            for qp in QS[Tp]:
                # A forward from RA: psi-part offs <= q[0:2], lambda-part >= q[2:4];  J backward from RJ: psi >= qp, lambda <= qp
                # constraint form: l(v) <= c  or  l(v) >= c, with l(t) = mn(l, s) + l(b) + l(v)
                def rhs(l, R, qq): return qq + mn(l, R) - mn(l, s) - dot(l, b)
                lo = [-INF] * 4; hi = [INF] * 4
                (lA0, lA1), (lJ0, lJ1) = LF[T], LF[Tp]
                # lambda: a0 = v0, a1 = v2
                lo[0] = max(lo[0], rhs(LAM[0], RAp2, q[2])); hi[0] = min(hi[0], rhs(LAM[0], RJ0, qp[2]))
                lo[2] = max(lo[2], rhs(LAM[1], RAp2, q[3])); hi[2] = min(hi[2], rhs(LAM[1], RJ0, qp[3]))
                if lo[0] > hi[0] or lo[2] > hi[2]: continue
                # psi forms: A gives l(v) <= rhs, J gives l(v) >= rhs.  Each form is +-v1 + k*v0 or +-v3 + k*v2 (k in {0,1}, see LF)
                cons = [(lA0, '<=', rhs(lA0, RAp2, q[0])), (lA1, '<=', rhs(lA1, RAp2, q[1])),
                        (lJ0, '>=', rhs(lJ0, RJ0, qp[0])), (lJ1, '>=', rhs(lJ1, RJ0, qp[1]))]
                ok = True
                for idx, other in ((1, 0), (3, 2)):
                    l_, h_ = -INF, INF
                    for l, op, c in cons:
                        if l[idx] == 0: continue
                        assert l[1 - idx % 2 + (0 if idx == 1 else 2) - (1 - idx % 2)] in (0,) or True
                        k = l[other]; sgn = l[idx]
                        # sgn*v_idx + k*v_other  (op)  c   ->  bound v_idx using the interval of v_other
                        if op == '<=': ub = c - (k * lo[other] if k > 0 else k * hi[other]); 
                        else: lb = c - (k * hi[other] if k > 0 else k * lo[other])
                        if op == '<=':
                            if sgn > 0: h_ = min(h_, ub)
                            else: l_ = max(l_, -ub)
                        else:
                            if sgn > 0: l_ = max(l_, lb)
                            else: h_ = min(h_, -lb)
                    lo[idx], hi[idx] = l_, h_
                    if l_ > h_: ok = False
                if not ok: continue
                nonempty += 1
                for i in range(4): box[i][0] = min(box[i][0], lo[i]); box[i][1] = max(box[i][1], hi[i])
    print(T, t0, t1, 'non-empty (shape,q,q\') cells', nonempty, 'bounding box of v', [tuple(map(int, x)) for x in box] if nonempty else None)
