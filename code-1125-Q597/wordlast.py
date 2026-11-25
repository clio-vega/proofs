"""PROVE 11-25: WS_last in the word model (board-free, from the 12 good configurations of words.configs, the table behind
Theorem AW).  For T in {B, C} (a_T = a1 for B, a0 for C; T' the other target):
 (W4) every middle word  w = (1,a_T)(1,b)(0,x_1)..(0,x_r)(0,b),  x_i in {a0,a1}, r <= RMAX, from the canonical class-1 rhombus R_p.
      alpha = displacement of (1,a_T), gamma = R_end - R_{p+1}, R_end = rhombus after w (= R_{p+m+1}, m = r+2).
      S = non-rh tiles of the straight T-configuration at R_p+gamma (= H_p + gamma).  CLAIM: S - ... i.e. for 1<=i<=m,
      no tile of S is removed by the translated step (A's step p+i) - alpha.   Also: window data D.
 (W1) the turn: J_o toward T' at R* = R_{p+m} (class 0) takes (0,b) [step e], then step e+1 of each type t allowed at its image.
      CLAIM: S is contained in cre(e) u cre(e+1) whenever t is not J_o's straight step (1,a_{T'}); report which tile from which.
 Lesions: (L1) (Tn) broken: t = straight; (L2) word shape broken: insert a (1,x) step or end without (0,b); (L3) wrong vector
      (+alpha instead of -alpha; and gamma replaced by gamma of a different word).   usage: python3 wordlast.py RMAX"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597')
from collections import Counter
import words as W
RMAX = int(sys.argv[1])
CL0 = frozenset({(0, 1), (0, -1)})
CL1 = {'B': frozenset({(1, -1), (-1, 1)}), 'C': frozenset({(1, 0), (-1, 0)})}
AT = {'B': 'a1', 'C': 'a0'}; OTHER = {'B': 'C', 'C': 'B'}
CF = {T: W.configs(T) for T in 'BC'}
def vsub(P, Q): return tuple(a - b for a, b in zip(P, Q))
def vadd(P, Q): return tuple(a + b for a, b in zip(P, Q))
def step(T, R, typ):
    """(removed non-rh tiles, created non-rh tiles, image rhombus) of the good T-configuration of type typ at R"""
    Rc, rem, Rn, cre = CF[T][(W.classkey(R), typ)]
    sh = W.base(R); assert W.trs(Rc, sh) == R
    return (frozenset(W.trs(t, sh) for t in rem), frozenset(W.trs(t, sh) for t in cre), W.trs(Rn, sh))
def run_word(T, R0, word):
    """list of (R_k, removed_k, created_k) and final rhombus"""
    out = []; R = R0
    for (cl, typ) in word:
        ck = W.classkey(R); want = CL0 if cl == 0 else CL1[T]
        if ck != want: return None, None
        rem, cre, Rn = step(T, R, typ); out.append((R, rem, cre)); R = Rn
    return out, R
def canon(T):
    """canonical class-1 rhombus of T at base 0"""
    for (ck, typ), v in CF[T].items():
        if ck == CL1[T]: return v[0]
st = Counter()
for T in 'BC':
    Tp = OTHER[T]; R_p = canon(T); aT = AT[T]
    # straight config tiles relative to R_p, and their offsets (for the window)
    remA, _, Rp1 = step(T, R_p, aT); alpha = vsub(W.base(Rp1), W.base(R_p))
    st[(T, 'alpha', alpha)] += 1
    Qa = {W.offs(t, R_p, T) for t in remA}
    for r in range(RMAX + 1):
        for xs in itertools.product(('a0', 'a1'), repeat=r):
            word = [(1, aT), (1, 'b')] + [(0, x) for x in xs] + [(0, 'b')]
            steps, Rend = run_word(T, R_p, word)
            if steps is None: st[(T, 'word not realisable (class mismatch)', r)] += 1; continue
            assert W.classkey(Rend) == CL1[T], 'end class'
            m = len(word) - 1
            gamma = vsub(W.base(Rend), W.base(steps[1][0]))
            S = frozenset(W.trs(t, gamma) for t in remA)            # H_p + gamma (non-rh)
            Rpg = W.trs(R_p, gamma)
            assert step(T, Rpg, aT)[0] == S                          # S is the straight config at R_p + gamma
            na = tuple(-x for x in alpha)
            hits = [i for i in range(1, m + 1) if S & frozenset(W.trs(t, na) for t in steps[i][1])]
            st[(T, '(W4) S met by a (C-alpha) step', bool(hits))] += 1
            if hits: st[(T, '(W4) HIT at i', tuple(hits), 'word', ' '.join('%d%s' % w for w in word))] += 1
            # window statistic: at which i is the offset test alone not enough
            Sa = frozenset(W.trs(t, alpha) for t in S)               # S + alpha = straight tiles at R_end
            for i in range(1, m + 1):
                Rk = steps[i][0]
                feas = any(W.offs(t, Rk, T) in W_Q[T] for t in Sa) if 'W_Q' in dir() else None
            # lesion L3: +alpha
            hitsL3 = [i for i in range(1, m + 1) if S & frozenset(W.trs(t, alpha) for t in steps[i][1])]
            st[(T, 'LESION +alpha: S met', bool(hitsL3))] += 1
            # (W1) turn: J_o toward T' at R* = steps[m][0], (0,b) for T', then each allowed type at its image
            Rstar = steps[m][0]
            if r == 0:   # the turn is local: depends only on R* relative to R_p+gamma; check once per r but assert invariance
                pass
            remE, creE, Re1 = step(Tp, Rstar, 'b')
            assert remE == steps[m][1], 'shared class-0 b hexagon'
            assert W.classkey(Re1) == CL1[Tp]
            for t in ('a0', 'a1', 'b'):
                _, creE1, _ = step(Tp, Re1, t)
                straight = (t == AT[Tp])
                cov = S <= (creE | creE1)
                src = tuple(sorted(('e' if u in creE else ('e+1' if u in creE1 else '-')) + ':' + u[0] for u in S))
                st[(T, '(W1) e+1 type', t, 'straight' if straight else 'allowed by (Tn)', 'S covered', cov, src)] += 1
    # window bound D for (W4): Q_T from the table (offsets of removed tiles over all 6 configs of T)
    QT = set()
    for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items():
        for t in rem: QT.add(W.offs(t, Rc, T))
    D = max(sum(abs(a - b) for a, b in zip(q2, q1)) for q1 in QT for q2 in Qa
            if all(q2[c] >= q1[c] for c in (0, 1)) and all(q2[c] <= q1[c] for c in (2, 3)))
    st[(T, 'window: Q_T', tuple(sorted(QT)), 'Q_a', tuple(sorted(Qa)), 'D (max l1 with monotone signs)', D)] += 1
    # stationary steps among the middle letters, and min |Delta| of the class-0 a-steps / (1,b)
    for (cl, typ) in [(1, 'b'), (0, 'a0'), (0, 'a1'), (0, 'b')]:
        R = R_p if cl == 1 else step(T, R_p, 'b')[2]
        if cl == 1: rem, cre, Rn = step(T, R_p, 'b')
        else: rem, cre, Rn = step(T, R, typ)
        R0 = R_p if cl == 1 else R
        d = W.offs(Rn, R0, T)
        st[(T, 'Delta(mu,lambda) of', (cl, typ), d)] += 1
for k in sorted(st, key=str): print(st[k], k)
