"""PROVE 11-23 c2, (W_a) last piece: for each J_o tile of steps c0..c0+2 not removed by A at p, p+1 (explicit, placed),
can A remove it at a step k >= p+2?  One-sided window (Lemma window): it would need q in Q_T with
om_T(t, R^A_{p+2}) <= q in the psi-part and >= q in the lambda-part.  Also all T-words of length <= 8 from R^A_{p+2}
(CONTROL that the window is not vacuous: the tiles removed by A at p+1 must pass the window from R^A_{p+1}).
usage: python3 wakeA3.py"""
import sys
sys.argv = ['wakeA.py', '2']
exec(open('wakeA.py').read().split('\nnr = lambda')[0])
nr = lambda S: {u for u in S if u[0] != 'rh'}
def winA(T, t, R):
    oA = offs(t, R, T)
    return [q for q in QS[T] if oA[0] <= q[0] and oA[1] <= q[1] and oA[2] >= q[2] and oA[3] >= q[3]]
for (T, t0, t1, off) in CASEA:
    Tp = OTHER[T]
    for t2 in ('a1', 'b') if T == 'B' else ('a0', 'b'):
        RJ0 = canon(Tp, t0[0]); H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1]); H2, RJ3, C2 = step(Tp, RJ2, t2)
        _, RAp, RAp1, RAp2, Hmid = place(T, t0, off)
        Hq = step(T, RAp1, STR[Tp])[0]
        left = (nr(H0) | nr(H1) | nr(H2)) - nr(Hmid)
        res = [(u[0], phi(u), 'J step', [i for i, H in enumerate((H0, H1, H2)) if u in H], 'window q', len(winA(T, u, RAp2))) for u in left]
        Sa = nr(SA(T, RAp2, 8))
        ctl = [len(winA(T, u, RAp1)) > 0 for u in nr(Hq)]
        print(T, t0, t1, t2, res, '| in words<=8 from R_p+2:', len(left & Sa), '| CONTROL p+1 tiles pass window from R_p+1:', ctl)
