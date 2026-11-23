"""PROVE 11-23 c2, Lemma grid: the termination test needs the packed vertex on nest T's 2-d grid, i.e.
f0_T(q) = f0_T(X_T) and f1_T(q) = f1_T(X_T), where f0_T is the Lambda_0 functional with kernel R u_T and f1_T the
Lambda_1 functional with kernel R w_T.  Frames (sheet tex Lemma fin): u_B=(-1,1;0,0), w_B=(0,0;1,-1), u_C=(1,0;0,0), w_C=(0,0;-1,0).
So f0_B = a0+b0, f1_B = a1+b1, f0_C = b0, f1_C = b1.  Table: Delta F (F = f(lexmin vertex)) over the 6 good configurations
of each target (translation invariant, so this is every step of every terminating journey by Lemma word).
LESION: the wrong target's functionals (f_C along B-journeys and vice versa) -- expect a sign change somewhere.
usage: python3 grid.py"""
import sys
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
F = {'B': (lambda P: P[0] + P[1], lambda P: P[2] + P[3]), 'C': (lambda P: P[1], lambda P: P[3])}
U = {'B': ((-1, 1, 0, 0), (0, 0, 1, -1)), 'C': ((1, 0, 0, 0), (0, 0, -1, 0))}
for T in 'BC':
    f0, f1 = F[T]; u, w = U[T]
    assert f0(u) == 0 and f1(w) == 0 and f0(w) == 0 and f1(u) == 0 and f0((1, 0, 0, 0)) != 0 or f0((0, 1, 0, 0)) != 0
    for G, tag in ((T, 'own'), (OTHER[T], 'LESION other target journeys')):
        rows = []
        for (ck, typ), (Rc, rem, Rn, cre) in sorted(CF[G].items(), key=str):
            cls = 0 if ck == k0 else 1
            d0 = f0(min(Rn[1])) - f0(min(Rc[1])); d1 = f1(min(Rn[1])) - f1(min(Rc[1]))
            rows.append(((cls, typ), int(d0), int(d1)))
        s0 = {r[1] for r in rows}; s1 = {r[2] for r in rows}
        print('functionals of', T, 'along', G, '-journeys (' + tag + '):', rows, '| Delta f0 values', sorted(s0), 'Delta f1 values', sorted(s1))
