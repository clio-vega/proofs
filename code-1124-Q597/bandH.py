"""Structural inputs at band vertices (all j in [c0+2,e], split records), in M0:
 (H1) every unit triangle of Delta at vbar_j has exactly one Lambda_k-lift in M0;
 (H2) across edge i (ccw from t_{i-1,i} to t_{i,i+1}) the height jump is k_i * n_i, k_i = #squares over E_i, with n_i
      depending only on (T, i), and n_{i+3} = -n_i, n_0 + n_2 = n_1;
 (H3) no rhombus of M0 lies over E_i;
 (Cl) closure k0-k3 = k4-k1 = k2-k5 (a consequence of H1, H2).   usage: python3 bandH.py SIZE"""
import sys
from collections import Counter, defaultdict
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/bandrule.py').read().split('st = Counter(); shapes')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
st = Counter(); NRM = defaultdict(set)
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTH[T]; k = ZKK[Tp]
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    lifts = defaultdict(list)
    for u in M0:
        if zisL(u, k): lifts[zpk(u, k)].append(u)
    SQc = Counter(zpk(u, k) for u in M0 if u[0] == "sq"); RHc = Counter(zpk(u, k) for u in M0 if u[0] == 'rh')
    for j in range(c0 + 2, e + 1):
        Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
        v = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1))); vb = zfree(v, k)
        H1 = all(len(lifts[tri_at(vb, i, 1)]) == 1 for i in range(6))
        K = []; H3 = True
        for i, d in enumerate(DIRS):
            E = frozenset((vb, (vb[0] + d[0], vb[1] + d[1]))); K.append(SQc[E]); H3 &= RHc[E] == 0
            if H1:
                h1 = zfix(min(lifts[tri_at(vb, i, -1)][0][1]), k); h2 = zfix(min(lifts[tri_at(vb, i, 1)][0][1]), k)
                jump = (h2[0] - h1[0], h2[1] - h1[1])
                if SQc[E]: NRM[(T, i)].add((jump[0] // SQc[E], jump[1] // SQc[E], jump[0] % SQc[E], jump[1] % SQc[E]))
                else: st[('zero-stack jump', jump)] += 1
        Cl = K[0] - K[3] == K[4] - K[1] == K[2] - K[5]
        st[(T, 'H1', H1, 'H3', H3, 'Cl', Cl)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
for kk in sorted(NRM): print('normal', kk, NRM[kk])
