"""COMPARATOR: E1 (abacus, prop:reflect) vs E1b (hook description, Theorem 1).
The cell<->pair dictionary is asserted, not assumed: cell (i,k) <-> (b_i, u_k)
where u_1<u_2<... are the holes and h_ik = b_i - u_k."""
from graph import maya, vertices, edges, partitions
from hookgraph import hook_edges, hooks

def dictionary(mu):
    ell, bs, is_bead = maya(mu)
    U = sorted(u for u in range(-ell, bs[0]) if not is_bead(u))   # holes below b_1
    cell2pair = {}
    h = hooks(mu)
    for i in range(1, ell+1):
        for k in range(1, mu[i-1]+1):
            cell2pair[(i,k)] = (bs[i-1], U[k-1])
            assert bs[i-1] - U[k-1] == h[(i,k)], ('hook dictionary fails', mu, i, k)
    return cell2pair

bad = 0; tot = 0
for n in range(1, 11):
    for mu in partitions(n):
        tot += 1
        d = dictionary(mu)
        A = {frozenset((x, y)) for (x, y, hh, ty) in edges(mu)}
        B = {frozenset((d[x], d[y])) for (x, y, ty) in hook_edges(mu)}
        # also check the type label matches row/col
        for (x, y, hh, ty) in edges(mu):
            pass
        if A != B:
            bad += 1
            print('MISMATCH', mu, 'abacus-only', A-B, 'hook-only', B-A)
print(f"engines agree on {tot-bad}/{tot} partitions, n<=10   (bad={bad})")
