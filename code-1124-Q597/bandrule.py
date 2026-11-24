"""The band induction step, as a finite configuration statement.  For band j in [c0+2, e], with T in {B,C}:
  (S) the two Lambda_k-triangles of Delta at edge s_T of vbar_j (s_B = 2, s_C = 3) both lift into H^o_{j-1} (M0 tiles)
      with equal heights;
  (I) the edge i_T (i_B = 3, i_C = 2) of vbar_j: if vbar_j = vbar_{j-1} + d, then E_{i_T}(vbar_j) = E_{i_T + 3}(vbar_{j-1})
      when d = DIRS[i_T+3], so the IH gives k = 0; at the base j = c0+2 both lifts lie in A's contact hexagon H_p.
Tabulate: offset o_j - o_{j-1}, type of step j-1, the projected step d, and whether (S), (I) hold, and whether the
offset/type determine (S).  Also controls: is the certificate's pair of lifts the same SHAPE (relative to R^o_{j-1})
for each (T, type_{j-1}, offset)?   usage: python3 bandrule.py SIZE"""
import sys
from collections import Counter, defaultdict
SIZE = int(sys.argv[1])
exec(open('/home/clio/projects/proofs/code-1124-Q597/vstar.py').read().split('st = Counter()')[0].replace("SIZE = int(sys.argv[1])", "SIZE = %d" % SIZE))
DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
SI = {'B': (2, 3), 'C': (3, 2)}
def tri_at(vb, i, s):  # triangle of Delta at vbar on side s=+1 (ccw) / -1 of edge i
    d = DIRS[i]; d2 = DIRS[(i + s) % 6]
    return frozenset((vb, (vb[0] + d[0], vb[1] + d[1]), (vb[0] + d2[0], vb[1] + d2[1])))
st = Counter(); shapes = defaultdict(set)
def rel(u, o): return (u[0], frozenset(tuple(a - b for a, b in zip(Q, o)) for Q in u[1]))
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if typed(Bq)[p][1] != 'b': continue
    T = 'BC'[L]; Tp = OTH[T]; k = ZKK[Tp]; TO = typed(Jo)
    M0 = set(P['M0']); W = M0 - set(apply(M0, Jo))
    rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
    over0 = {zpk(u, k): u for u in M0 if zisL(u, k)}
    def V(j):
        Rj = next(u for u in Jo[j][1] if u[0] == 'rh'); o = min(Rj[1])
        return o, (o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1))))
    sE, iE = SI[T]
    for j in range(c0 + 2, e + 1):
        oj, vj = V(j); op, vp = V(j - 1); vb = zfree(vj, k); vbp = zfree(vp, k)
        d = (vb[0] - vbp[0], vb[1] - vbp[1])
        H = Jo[j - 1][1] & M0; overH = {zpk(u, k): u for u in H if zisL(u, k)}
        a, b = tri_at(vb, sE, 1), tri_at(vb, sE, -1)
        S = a in overH and b in overH and zfix(min(overH[a][1]), k) == zfix(min(overH[b][1]), k)
        if S: shapes[(T, TO[j - 1], tuple(x - y for x, y in zip(oj, op)))].add((rel(overH[a], op), rel(overH[b], op)))
        if j == c0 + 2:
            Hp = A[p][1] & M0; overP = {zpk(u, k): u for u in Hp if zisL(u, k)}
            a2, b2 = tri_at(vb, iE, 1), tri_at(vb, iE, -1)
            I = a2 in overP and b2 in overP and zfix(min(overP[a2][1]), k) == zfix(min(overP[b2][1]), k)
            how = 'base:Hp'
        elif d == (0, 0): I = True; how = 'same vertex'
        else: I = (d == DIRS[(iE + 3) % 6]); how = 'moved %s' % (d,)
        st[(T, 'type j-1', TO[j - 1], 'off', tuple(x - y for x, y in zip(oj, op)), 'd', d, '(S)', S, '(I)', I, how)] += 1
for kk in sorted(st, key=str): print(st[kk], kk)
print('distinct certificate shapes per (T, type_{j-1}, offset):')
for kk in sorted(shapes, key=str): print(' ', kk, len(shapes[kk]))
# --- shapes, printed (relative to o_{j-1}) ---
for kk in sorted(shapes, key=str):
    for sh in shapes[kk]: print('SHAPE', kk, [ (x[0], sorted(x[1])) for x in sh])
