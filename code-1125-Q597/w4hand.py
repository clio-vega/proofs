"""(W4) by hand-size checks.  Offsets om(t,R) = (psi(t)-mu(R), lambda(t)-lambda(R)) [words.offs].
 (a) Q_T, Q_a (offsets of the straight config's removed non-rh tiles at R_p), Delta of (1,b), (0,a0), (0,a1), (0,b).
 (b) CONE: is there q_a in Q_a, q in Q_T, delta in {Delta(1,b), Delta(0,a0), Delta(0,a1)} with q - q_a - delta in the cone
     K = {mu-part <= 0, lambda-part >= 0}?  (If none: no step before the last of the middle removes a tile of S + alpha.)
 (c) LOCAL: rem(T-(0,b) at R*) meets rem(straight T-config at the image of that (0,b))?  (must be empty)
 Lesion for (b): drop delta (allow Sigma = 0) -> must become feasible (that is the i = m case, handled by (c))."""
from wmodel import *
def cone(d): return d[0] <= 0 and d[1] <= 0 and d[2] >= 0 and d[3] >= 0
for T in 'BC':
    R_p = canon(T); aT = AT[T]
    remA, _, Rp1 = step(T, R_p, aT)
    QT = set()
    for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items():
        for t in rem: QT.add(W.offs(t, Rc, T))
    Qa = sorted({W.offs(t, R_p, T) for t in remA})
    R1 = Rp1; _, _, R0c = step(T, R1, 'b')         # (1,b) from R_{p+1} to a class-0 rhombus
    D = {}
    D['1b'] = W.offs(R0c, R1, T)
    for typ in ('a0', 'a1', 'b'):
        D['0' + typ] = W.offs(step(T, R0c, typ)[2], R0c, T)
    print(T, 'Q_T', sorted(QT)); print(T, 'Q_a', Qa); print(T, 'Delta', D)
    feas = [(qa, q, k) for qa in Qa for q in QT for k in ('1b', '0a0', '0a1')
            if cone(tuple(a - b - c for a, b, c in zip(q, qa, D[k])))]
    print(T, '(b) CONE feasible triples:', feas)
    les = [(qa, q) for qa in Qa for q in QT if cone(tuple(a - b for a, b in zip(q, qa)))]
    print(T, '    lesion (no delta): feasible pairs', les)
    # (c) local
    remE, _, Rend = step(T, R0c, 'b')
    remS, _, _ = step(T, Rend, aT)
    print(T, '(c) LOCAL: (0,b) hexagon meets straight config at its image:', bool(remE & remS))
    # (c) lesion: straight config at the image of (0,a0)/(0,a1)? not class 1 -> use the T-(0,b) image's other a-type
    other = 'a0' if aT == 'a1' else 'a1'
    print(T, '    lesion: (0,b) hexagon meets the NON-straight class-1 a-config / b-config at its image:',
          bool(remE & step(T, Rend, other)[0]), bool(remE & step(T, Rend, 'b')[0]))
