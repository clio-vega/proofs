"""(Ea) why: for R in (R_p, R_{p+1}, R'_{p+1}, R_{p+2}) in case a: rhombus edge set (Lambda_0 edge, Lambda_1 edge) vs nest T's
frame (u_T, w_T); and for the lexmin vertex q: is (q - X_T)'s Lambda_0 block on the line R u_T, its Lambda_1 block on R w_T?
Same statistics at every rhombus of every journey (CONTROL / context), and at final rhombi.   usage: python3 eaB.py SIZE"""
import sys
exec(open('eaA.py').read().split("\naT = ('a1', 'a0')")[0])
def blocks(info, R):
    S = sorted(R if not isinstance(R, tuple) else R[1]); q = S[0]
    E = {tuple(b - a for a, b in zip(S[0], P)) for P in S[1:]}
    u, w = info['u'], info['w']
    e0 = [e for e in E if not any(e[2:])]; e1 = [e for e in E if not any(e[:2])]
    shape = (e0 and (tuple(e0[0]) in (u, tuple(-x for x in u)))) and (e1 and (tuple(e1[0]) in (w, tuple(-x for x in w))))
    D = tuple(a - b for a, b in zip(q, info['X']))
    on0 = coef((D[0], D[1], 0, 0), u) is not None if any(D[:2]) else True
    on1 = coef((0, 0, D[2], D[3]), w) is not None if any(D[2:]) else True
    return ('shape', bool(shape), 'L0 on line', on0, 'L1 on line', on1)
aT = ('a1', 'a0'); st = Counter()
for P in PAIRS:
    N0 = P['nests']
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    TT = 1 if L == 0 else 2
    if typed(Bq)[p][1] == 'b':
        st[('case b ctx', 'T', 'BC'[L], 'R_p', blocks(N0[TT], rh(A[p][1])))] += 1
        continue
    for nm, R in (('R_p', rh(A[p][1])), ('R_p+1', rh(A[p][2])), ("R'_p+1", rh(Bq[p][2])), ('R_p+2', rh(A[p+1][2]))):
        st[('T', 'BC'[L], nm, blocks(N0[TT], R))] += 1
    st[('T', 'BC'[L], 'final', blocks(N0[TT], rh(A[-1][2])))] += 1
    for i in range(len(A)):
        st[('T', 'BC'[L], 'all rhombi of A', 'class', typed(A)[i][0], blocks(N0[TT], rh(A[i][1])))] += 1
for k in sorted(st, key=str): print(st[k], k)
