"""Compare E2 (matrix, rim hooks) against E1 (abacus, prop:reflect).
Column gamma of M <-> vertex (b,b') of G: gamma = B(mu) with bead b moved to hole b'.
The dictionary is COMPUTED from gamma, not assumed."""
from graph import maya, vertices, edges, partitions
from matrix_engine import two_term_relations

def col2pair(mu, gam):
    ell, bs, is_bead = maya(mu)
    Bmu = set(bs) | {-k for k in range(ell+1, ell+sum(mu)+5)}
    ellg = len(gam)
    # use the SAME ambient normalisation: pad gamma to the same number of rows via mu_i=0 -> -i
    N = ell + sum(mu) + 5
    Bm = {mu[i]-(i+1) if i < ell else -(i+1) for i in range(N)}
    Bg = {(gam[i] if i < ellg else 0)-(i+1) for i in range(N)}
    out = Bm - Bg; inn = Bg - Bm
    assert len(out) == 1 and len(inn) == 1, (mu, gam, out, inn)
    return (out.pop(), inn.pop())          # (bead b, hole b')

bad = 0; tot = 0; nvtx_mismatch = 0
for n in range(1, 9):
    for mu in partitions(n):
        tot += 1
        cols, rels = two_term_relations(mu)
        V = set(vertices(mu))
        P = {g: col2pair(mu, g) for g in cols}
        assert set(P.values()) == V, (mu, 'vertex set differs')
        A = {frozenset((x, y)): h for (x, y, h, ty) in edges(mu)}
        B = {}
        for (g1, g2, h) in rels:
            k = frozenset((P[g1], P[g2]))
            B.setdefault(k, set()).add(h)
        # edge sets equal?
        if set(A) != set(B):
            bad += 1; print('EDGE MISMATCH', mu, 'reflect-only', set(A)-set(B), 'matrix-only', set(B)-set(A))
            continue
        # labels equal?  a matrix row may present the relation in either direction (+h or -h)
        for k in A:
            hs = B[k]
            if not (A[k] in hs or -A[k] in hs):
                bad += 1; print('LABEL MISMATCH', mu, k, 'reflect', A[k], 'matrix', hs)
print(f"E1 (prop:reflect) vs E2 (matrix two-term rows): agree on {tot-bad}/{tot} partitions, n<=8")
