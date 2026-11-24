"""Draw J_o's board at step e (before the split step) near R*, for one case-b record of target T, with: HJ (split hexagon,
grey), R* (black), J_o's straight config at R^o_{e+1} (red outline), A's actual next hexagon at R# (blue outline),
J_o's actual step-e+1 hexagon (green outline), universal blockers (orange fill).  usage: python3 draw.py SIZE T idx"""
import sys, os, pickle
SIZE, TT, IDX = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
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
    T = 'BC'[L]; Tp = OTHER[T]
    if T != TT or typed(Bq)[p][1] != 'b': continue
    n += 1
    if n == IDX: break
M0 = set(P['M0']); Mx = set(apply(M0, Jo)); W = M0 - Mx
rem = {}
for j, s in enumerate(Jo):
    for u in s[1]: rem.setdefault(u, j)
c0 = min(rem[u] for u in A[p][1] if u in W); e = c0 + len(mA) - 1
bd = set(M0)
for s in Jo[:e]: bd = (bd - s[1]) | s[2]
HJ = Jo[e][1]; Rs = next(u for u in HJ if u[0] == 'rh'); RJ = next(u for u in Jo[e][2] if u[0] == 'rh')
S, _, _ = step(Tp, nz(RJ), STRT[Tp])
bd1 = (bd - Jo[e][1]) | Jo[e][2]
nb = {nz(u): u for u in bd1}; blk = {nb[b] for (s, b) in conflict(S, set(nb))}
c = phi(min(Rs[1])); fig, ax = plt.subplots(figsize=(9, 9))
def draw(tile, **kw):
    zs = poly(nz(tile)); ax.fill([z.real for z in zs], [z.imag for z in zs], **kw)
for u in bd:
    zs = poly(nz(u))
    if min(abs(z - c) for z in zs) > 4.5: continue
    fc = {'tri': '#eeeeee', 'sq': '#dde8ff', 'rh': '#ffe8cc'}[u[0]]
    if u in HJ: fc = '#aaaaaa'
    if u in blk: fc = 'orange'
    draw(u, facecolor=fc, edgecolor='k', lw=0.5)
    m = sum(zs)/len(zs)
    tag = ('W' if u in W else ('X' if u in M0 else 'j%d' % next((j for j, s in enumerate(Jo) if u in s[2]), -1)))
    ax.text(m.real, m.imag, tag, fontsize=6, ha='center')
draw(Rs, facecolor='black', alpha=0.5)
for u in S: draw(u, facecolor='none', edgecolor='red', lw=2.5, ls='--')
if e + 1 < len(Jo):
    for u in Jo[e+1][1]: draw(u, facecolor='none', edgecolor='green', lw=2)
k = p + len(mA)
if k < len(A):
    for u in A[k][1]: draw(u, facecolor='none', edgecolor='blue', lw=2)
ax.set_aspect('equal'); ax.set_title('T=%s rec %d  m=%d  Jo e+1 type %s  A ends at R#: %s' % (T, IDX, len(mA)-1, typed(Jo)[e+1], k >= len(A)))
plt.savefig('/home/clio/projects/proofs/code-1124-Q597/draw-%s-%d.png' % (T, IDX), dpi=80)
print('saved')
