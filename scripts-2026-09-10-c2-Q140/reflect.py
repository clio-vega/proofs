"""Two-term rows, via the abacus.  A row lambda != mu with exactly two nonzero
entries gives a relation  t^h w(col1) + t^h' w(col2) = 0.

DERIVED (to be checked against the matrix):
 (II) same bead b, holes b'<e<b, c:=b'+e-b in B and c<b':
        w(b,e) = -t^{1+|B cap (b',e)|} w(b,b')
 (I)  same hole b', beads c<b, e:=b+c-b' not in B and e>b:
        w(c,b') = -t^{|B cap (c,b)|} w(b,b')
"""
from rim import *

def maya(mu, N):
    return set(beta_set(mu, N))

def cols_abacus(mu, N):
    B = maya(mu, N)
    out = {}
    for b in sorted(B):
        for bp in range(0, b):
            if bp not in B:
                g = from_beta((B - {b}) | {bp})
                out[(b, bp)] = g
    return B, out

def derived_relations(mu, N):
    B, cols = cols_abacus(mu, N)
    rels = []
    for (b, bp) in cols:
        for e in range(bp+1, b):
            if e in B: continue
            c = bp + e - b
            if c in B and c < bp:                      # type (II)
                h = 1 + len([x for x in B if bp < x < e])
                rels.append((("II"), (b,bp), (b,e), h))   # w(b,e) = -t^h w(b,bp)
    for (b, bp) in cols:
        for c in sorted(B):
            if c >= b or c <= bp: continue
            e = b + c - bp
            if e not in B and e > b:                    # type (I)
                h = len([x for x in B if c < x < b])
                rels.append((("I"), (b,bp), (c,bp), h))   # w(c,bp) = -t^h w(b,bp)
    return B, cols, rels

# ---- calibrate the derived relations against the matrix, then test sufficiency
print(" n | mu                  | #2-term rows | rank(2-term only) | rank(derived) | n")
for n in range(4, 9):
    for mu in partitions(n):
        N = n + len(mu) + 2
        rows, cols, M, b = local_system(mu)
        # empirical two-term rows
        rel_rows = [i for i,l in enumerate(rows)
                    if l != tuple(mu) and sum(1 for j in range(M.cols) if M[i,j]!=0) == 2]
        R = M[rel_rows, :] if rel_rows else sp.zeros(0, M.cols)
        r_emp = R.rank() if rel_rows else 0
        # derived relations, as a matrix in the same column order
        B, cabac, rels = derived_relations(mu, N)
        colpos = {tuple(g): j for j, g in enumerate(cols)}
        D = []
        for kind, c1, c2, h in rels:
            g1, g2 = cabac[c1], cabac[c2]
            row = [0]*M.cols
            row[colpos[g1]] = t**h; row[colpos[g2]] = 1
            D.append(row)
        D = sp.Matrix(D) if D else sp.zeros(0, M.cols)
        r_der = D.rank() if D.rows else 0
        print(f"{n:2d} | {str(mu):19s} | {len(rel_rows):3d} | {r_emp:2d} {'FULL' if r_emp==n else '    '} "
              f"| {r_der:2d} {'FULL' if r_der==n else '    '} | {n}")
