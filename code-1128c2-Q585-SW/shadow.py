"""Shadow-then-merge test in puzzle coordinates.  K,H = x's solo B/C gash paths (positions incl. start), G = y->C after K,
Gp = y->B after H.  D = start(x) - start(y).  Classify:
 B-side: Gp = (K - D)[0..t] ++ K[t..]   (shadow t+1 positions, one connecting step, then x's path from index t)
 C-side: G  = (H - D)[0..t-1] ++ [one free position] ++ H[t+1..]  or other.
usage: python3 shadow.py NSEED0 NSEED1"""
import sys
src = open('/home/clio/projects/proofs/code-1128c2-Q585-SW/gpath.py').read().split("N0, N1 =")[0]
exec(src)
def sh(P, D): return [[(a - D[0], b - D[1]) for (a, b) in p] for p in P]
def classifyB(K, Gp, D):
    KD = sh(K, D)
    for t in range(len(K)):
        if Gp[:t+1] == KD[:t+1] and Gp[t+1:] == K[t:]: return ('shadow-then-K', t)
    return ('other',)
def classifyC(H, G, D):
    HD = sh(H, D)
    for t in range(len(H)):
        if G[:t] == HD[:t] and G[t+1:] == H[t+1:] and len(G) == len(H): return ('shadow-free-H', t, 'free==HD' if G[t] == HD[t] else ('free==H' if G[t] == H[t] else 'free new'))
    s = 0
    while s < min(len(G), len(H)) and G[s] == HD[s]: s += 1
    return ('other', 'shadowlen', s, 'lenG-lenH', len(G) - len(H), 'end equal', G[-1] == H[-1])
st = Counter(); S0, S1 = int(sys.argv[1]), int(sys.argv[2])
for k, (key, b0) in enumerate(census()):
    if k < S0 or k >= S1: continue
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            b1 = copy.deepcopy(b); K, _ = path(b1, x, dirs['B']); G, _ = path(b1, y, dirs['C'])
            b2 = copy.deepcopy(b); H, _ = path(b2, x, dirs['C']); Gp, _ = path(b2, y, dirs['B'])
            D = (K[0][0][0] - G[0][0][0], K[0][0][1] - G[0][0][1])
            st['n'] += 1; st[('D', D == (H[0][0][0]-Gp[0][0][0], H[0][0][1]-Gp[0][0][1]))] += 1
            cb = classifyB(K, Gp, D); cc = classifyC(H, G, D)
            st[('B',) + cb[:1]] += 1; st[('C',) + cc[:1] + cc[2:3]] += 1
            if cc[0] == 'other': st[('C other',) + cc[1:]] += 1
            st[('B t', cb[1] if len(cb) > 1 else None, 'C t', cc[1] if cc[0] != 'other' else None)] += 1
            st[('pre(K,H)', next((i for i in range(min(len(K), len(H))) if K[i] != H[i]), None), 'Bt', cb[1] if len(cb) > 1 else None)] += 1
for k in sorted(st, key=str): print(st[k], k)
