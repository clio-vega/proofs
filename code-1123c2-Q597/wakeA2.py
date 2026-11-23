"""PROVE 11-23 c2, (W_a) completion.
(1) FR refined: g(R^J_j) for j = c0, c0+1, c0+2, c0+3 (J's three placed steps), and the non-rh tiles of J's steps c0..c0+2
    not removed by A at p or p+1 (they are the only J tiles that the g-bound does not exclude when the gap is 0).
(2) the past window at a larger box (adequacy), and LESION: A's reference R_p (steps p, p+1 included) and J's reference
    R^J_{c0+3} (steps c0..c0+2 counted as past): feasible tiles must now exist, and the words must meet.
usage: python3 wakeA2.py [box]"""
import sys
box = sys.argv[1] if len(sys.argv) > 1 else '6'
sys.argv = ['wakeA.py', box]
exec(open('wakeA.py').read().split('\nnr = lambda')[0])
nr = lambda S: {u for u in S if u[0] != 'rh'}
for (T, t0, t1, off) in CASEA:
    Tp = OTHER[T]
    for t2 in ('a0', 'a1', 'b'):
        RJ0 = canon(Tp, t0[0]); H0, RJ1, C0 = step(Tp, RJ0, t0[1]); H1, RJ2, C1 = step(Tp, RJ1, t1[1])
        try: H2, RJ3, C2 = step(Tp, RJ2, t2)
        except KeyError: print(T, t0, t1, t2, 'not placeable'); continue
        _, RAp, RAp1, RAp2, Hmid = place(T, t0, off)
        Jrem = nr(H0) | nr(H1) | nr(H2)
        if len(nr(step(T, RAp, STR[T])[0]) & Jrem) != 3: print(T, t0, t1, t2, '|S| != 3, not a survivor'); continue
        left = sorted((phi(u), u[0]) for u in Jrem - nr(Hmid))
        gs = [g(R) for R in (RJ0, RJ1, RJ2, RJ3)]
        # A's tiles at k >= p+2: T=B phi >= g(R_p+2); T=C phi <= g(R_p+2).  J tiles at j: T=B (J to C) phi <= g(R^J_j); T=C phi >= g(R^J_j)
        if T == 'B': sep = [gs[j] < g(RAp2) for j in range(4)]; leftok = all(f < g(RAp2) for f, _ in left)
        else: sep = [gs[j] > g(RAp2) for j in range(4)]; leftok = all(f > g(RAp2) for f, _ in left)
        print(T, t0, t1, 'c0+2', t2, 'g(R^J_c0..c0+3)', gs, 'g(R_p+2)', g(RAp2), 'separated at j:', sep,
              '| J tiles of c0..c0+2 not removed at p,p+1 (phi,kind):', left, 'excluded by phi:', leftok)
        if t2 == 'b':
            DA, DJ, nf, ob = window(T, Tp, RAp2, RJ0)
            print('   past window box', BOX, 'feasible', nf, 'on bd', ob)
            DA, DJ, nf, ob = window(T, Tp, RAp, RJ3)
            print('   LESION window (A ref R_p, J ref R^J_c0+3): feasible', nf, 'on bd', ob, 'D_A', DA, 'D_J', DJ)
