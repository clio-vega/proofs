"""
Q149 engine: ribbon operators with FREE *SHAPE* weights.

The weight of the ribbon created by the bead move b -> b+e is an opaque symbol
indexed by the OCCUPANCY WORD u in {0,1}^{e-1} of the open interval (b, b+e),
equivalently by the composition alpha |= e of its row lengths.

Nothing here hard-codes (-1)^ht, t^ht, or height-only dependence.

Two independent ribbon enumerations, cross-checked:
  * ribbons_young : brute force over partitions, test border strip, read off the
                    row-length composition of the skew shape directly
  * ribbons_maya  : bead move, read off the occupancy word of (b, b+e)
"""
from itertools import combinations
import sympy as sp

# ---------------- partitions ----------------
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,)+rest

def contains(mu, lam):
    if len(lam) > len(mu): return False
    return all(mu[i] >= lam[i] for i in range(len(lam)))

def cells(lam):
    return {(i,j) for i,p in enumerate(lam) for j in range(p)}

# ---------------- word <-> composition ----------------
def word_to_comp(u):
    """u in {0,1}^{e-1}  ->  the row-length composition of the ribbon,
    read TOP ROW FIRST (the Young-diagram convention).

    The gaps between consecutive occupied sites of {b} u {b+i : u_i=1} u {b+e}
    are the row lengths read from the BOTTOM row up, because in beta-coordinates
    a larger bead sits in an earlier (higher) row.  Hence the reversal.
    Verified against ribbons_young by _crosscheck()."""
    occ = [0] + [i+1 for i,x in enumerate(u) if x] + [len(u)+1]
    return tuple(occ[i+1]-occ[i] for i in range(len(occ)-1))[::-1]

def comp_to_word(alpha):
    e = sum(alpha); u = [0]*(e-1); s = 0
    for a in alpha[::-1][:-1]:
        s += a; u[s-1] = 1
    return tuple(u)

def all_words(k):
    return [tuple((i>>j)&1 for j in range(k)) for i in range(2**k)]

# ---------------- implementation 1: Young diagram ----------------
def ribbons_young(lam, e):
    """list of (mu, rowcomp) where rowcomp is the composition of e given by the
    row lengths of the skew ribbon mu/lam, read TOP row first."""
    out = []; n = sum(lam)+e; L = cells(lam)
    for mu in partitions(n):
        if not contains(mu, lam): continue
        S = cells(mu)-L
        if len(S) != e: continue
        if any((i,j) in S and (i+1,j) in S and (i,j+1) in S and (i+1,j+1) in S for (i,j) in S):
            continue
        start = next(iter(S)); seen = {start}; stack=[start]
        while stack:
            (i,j) = stack.pop()
            for (p,q) in ((i+1,j),(i-1,j),(i,j+1),(i,j-1)):
                if (p,q) in S and (p,q) not in seen:
                    seen.add((p,q)); stack.append((p,q))
        if len(seen) != e: continue
        rows = sorted({i for (i,j) in S})
        comp = tuple(len([1 for (i,j) in S if i==r]) for r in rows)
        out.append((mu, comp))
    return out

# ---------------- implementation 2: Maya / beta-set ----------------
def maya(lam, L):
    ext = list(lam)+[0]*(L+2)
    return frozenset(ext[j-1]-j for j in range(1,L+2) if ext[j-1]-j >= -L)

def unmaya(M, L):
    beads = sorted(M, reverse=True); lam=[]
    for j,b in enumerate(beads, start=1):
        p = b+j
        if p > 0: lam.append(p)
        else: break
    return tuple(lam)

def ribbons_maya(lam, e, L=None):
    """list of (mu, word) with word = occupancy of (b, b+e)."""
    if L is None: L = sum(lam)+e+4
    M = maya(lam, L); out=[]
    for b in sorted(M):
        if (b+e) in M: continue
        M2 = (M-{b})|{b+e}
        u = tuple(1 if (b+i) in M else 0 for i in range(1, e))
        out.append((unmaya(M2,L), u))
    return out

def _crosscheck(maxn=8, maxe=5, verbose=True):
    """Young rowcomp  ==  word_to_comp(maya word), for every lam, e."""
    checked = 0
    for n in range(0, maxn+1):
        for lam in partitions(n):
            for e in range(1, maxe+1):
                a = sorted((mu, comp) for mu, comp in ribbons_young(lam, e))
                b = sorted((mu, word_to_comp(u)) for mu, u in ribbons_maya(lam, e))
                assert a == b, (lam, e, a, b)
                checked += len(a)
    if verbose: print(f"cross-check young(rowcomp) vs maya(word->comp): OK, {checked} ribbons")
    return True

def _crosscheck_word_comp(maxe=7):
    for e in range(1, maxe+1):
        ws = all_words(e-1)
        cs = [word_to_comp(u) for u in ws]
        assert len(set(cs)) == 2**(e-1), e
        assert all(sum(c)==e for c in cs)
        assert all(comp_to_word(word_to_comp(u))==u for u in ws)
    print(f"word<->composition bijection: OK for e<= {maxe}")
    return True

if __name__ == "__main__":
    _crosscheck_word_comp()
    _crosscheck()
