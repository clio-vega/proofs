"""Verify the bead-move <-> rim-hook dictionary at the level of Young diagrams,
and identify which orientation gives Postnikov's h_r (horizontal) and which e_r
(vertical)."""
from itertools import combinations

def word_to_partition(S, n, k):
    """S = set of occupied sites (0-indexed) in 0..n-1, |S|=k.
    beta sorted decreasing = lambda_j + k - j, j=1..k."""
    beta = sorted(S, reverse=True)
    lam = [beta[j] - (k - 1 - j) for j in range(k)]
    return tuple(x for x in lam if x > 0)

def diagram(lam):
    return set((r+1, c+1) for r, p in enumerate(lam) for c in range(p))

def rim_hook_data(lam, mu):
    """lam subset mu, mu/lam should be a rim hook.  Return (size, rows, cols,
    is_horizontal_strip, is_vertical_strip, connected, no2x2)."""
    D = diagram(mu) - diagram(lam)
    rows = set(r for r, c in D); cols = set(c for r, c in D)
    no2x2 = all(not ({(r,c),(r+1,c),(r,c+1),(r+1,c+1)} <= D) for r,c in D)
    # connectivity (edge-adjacent)
    if D:
        seen={next(iter(D))}; stack=list(seen)
        while stack:
            r,c=stack.pop()
            for nb in ((r+1,c),(r-1,c),(r,c+1),(r,c-1)):
                if nb in D and nb not in seen:
                    seen.add(nb); stack.append(nb)
        conn = (seen==D)
    else:
        conn=True
    horiz = len(cols)==len(D)   # at most one box per column
    vert  = len(rows)==len(D)   # at most one box per row
    return len(D), len(rows), len(cols), horiz, vert, conn, no2x2

def check(n, kmax=None):
    bad=[]
    for k in range(0, n+1):
        for Sc in combinations(range(n), k):
            S=set(Sc)
            lam = word_to_partition(S, n, k)
            for j in S:
                for e in range(1, n):
                    tgt = j+e
                    if tgt >= n or tgt in S:      # NON-CYCLIC moves only
                        continue
                    T = (S-{j})|{tgt}
                    mu = word_to_partition(T, n, k)
                    h = sum(1 for s in range(1,e) if j+s in S)
                    size, rows, cols, horiz, vert, conn, no2x2 = rim_hook_data(lam, mu)
                    ok = (size==e) and conn and no2x2 and (rows-1==h) and (cols-1==e-1-h)
                    if not ok:
                        bad.append((n,k,tuple(sorted(S)),j,e,lam,mu,h,size,rows,cols))
    return bad

for n in range(2,10):
    b=check(n)
    print(f'n={n}: violations of [bead move j->j+e over h beads == rim hook of size e, leg length (rows-1) = h, arm = e-1-h]:', len(b))
    for x in b[:3]: print('   ',x)
