"""Draw M0 near the whole shadow of one split record with m >= MM: M0 tiles (W shaded pink, Mx white), J_o's rhombi
R^o_j for j in [c0, e+1] (labelled), A's rhombi R_p..R_{p+m+1} (blue labels), H_p outlined blue, the two universal
triangles Z1, Z2 (orange), the straight configuration S at R^o_{e+1} (red dashed).  usage: python3 draw2.py SIZE T MM idx"""
import sys, os, pickle
SIZE, TT, MM, IDX = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
os.chdir('/home/clio/projects/proofs/code-1121-Q597')
exec(open("diverge.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import conflict, poly, phi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
STRT = {'B': 'a1', 'C': 'a0'}; aT = ('a1', 'a0')
def nz(t): return (t[0], frozenset(tuple(Fr(x) for x in P) for P in t[1]))
n = -1
for P in PAIRS:
    if P['kind'] != 'BC': continue
    T_ = [typed(A)[middles(A, Bq)[2]] for (t, A, Bq, Jo) in P['roles']]
    L = [i for i in (0, 1) if T_[i] == (1, aT[i])][0]
    t, A, Bq, Jo = P['roles'][L]; mA, mB, p, q = middles(A, Bq)
    if 'BC'[L] != TT or typed(Bq)[p][1] != 'b' or len(mA) - 1 < MM: continue
    n += 1
    if n == IDX: break
T = TT; Tp = OTHER[T]; TO = typed(Jo)
M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
rem = {}
for j, s in enumerate(Jo):
    for u in s[1]: rem.setdefault(u, j)
c0 = min(rem[u] for u in A[p][1] if u in W); m = len(mA) - 1; e = c0 + m
bd = set(M0)
for s in Jo[:e+1]: bd = (bd - s[1]) | s[2]
RJ = next(u for u in Jo[e][2] if u[0] == 'rh'); S, _, _ = step(Tp, nz(RJ), STRT[Tp])
nb = {nz(u): u for u in bd}; Z = {nb[b] for (_, b) in conflict(S, set(nb)) if nb[b] in M0 and b[0] == 'tri'}
fig, ax = plt.subplots(figsize=(10, 10))
def draw(tile, **kw):
    zs = poly(nz(tile)); ax.fill([z.real for z in zs], [z.imag for z in zs], **kw); return sum(zs)/len(zs)
for u in M0:
    fc = '#ffd0d0' if u in W else ('#ffe8cc' if u[0] == 'rh' else 'white')
    if u in Z: fc = 'orange'
    c = draw(u, facecolor=fc, edgecolor='grey', lw=0.5)
    if u in W: ax.text(c.real, c.imag, str(rem[u] - c0), fontsize=6, ha='center', color='red')
for j in range(c0, min(e + 2, len(Jo))):
    R = next(u for u in Jo[j][1] if u[0] == 'rh'); c = draw(R, facecolor='none', edgecolor='green', lw=1.5)
    ax.text(c.real, c.imag + 0.1, 'J%d' % (j - c0), fontsize=8, color='green')
for k in range(p, min(p + m + 2, len(A))):
    R = next(u for u in A[k][1] if u[0] == 'rh'); c = draw(R, facecolor='none', edgecolor='blue', lw=1.5)
    ax.text(c.real, c.imag - 0.15, 'A%d' % (k - p), fontsize=8, color='blue')
for u in A[p][1]: draw(u, facecolor='none', edgecolor='blue', lw=2.5)
for u in S: draw(u, facecolor='none', edgecolor='red', lw=2, ls='--')
ax.set_aspect('equal'); ax.set_title('T=%s m=%d  J_o word c0..e+1: %s' % (T, m, ' '.join('%d%s' % w for w in TO[c0:e+2])))
plt.savefig('/home/clio/projects/proofs/code-1124-Q597/draw2-%s-%d-%d.png' % (T, MM, IDX), dpi=75)
print('ok')
