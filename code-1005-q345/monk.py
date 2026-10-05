"""Monk transfer matrices M_p, and the H2 ground-truth check against
the divided-difference Schubert instrument (a genuinely different mechanism)."""
from itertools import permutations
from schub import length, rmul_t, schubert_table, pmul, extract_coeff

def monk_matrix(n, p):
    """M_p[u][v] = 1 iff v = u t_{ab} is a Bruhat cover with a <= p < b."""
    perms = list(permutations(range(1, n+1)))
    idx = {w: i for i, w in enumerate(perms)}
    M = [[0]*len(perms) for _ in perms]
    for u in perms:
        lu = length(u)
        for a in range(1, p+1):
            for b in range(p+1, n+1):
                v = rmul_t(u, a, b)
                if length(v) == lu+1:
                    M[idx[u]][idx[v]] += 1
    return M, perms, idx

def check_monk(n):
    """Verify Monk's rule against the Schubert-polynomial instrument."""
    N = n+1
    tab = schubert_table(n, N)
    perms = list(permutations(range(1, n+1)))
    sp = [tuple(list(range(1, p))+[p+1, p]+list(range(p+2, n+1))) for p in range(1, n)]
    bad = tested = 0
    for p in range(1, n):
        M, _, idx = monk_matrix(n, p)
        for u in perms:
            P = pmul(tab[sp[p-1]], tab[u])
            for v in perms:
                if length(v) != length(u)+1:
                    continue
                tested += 1
                truth = extract_coeff(P, v, N)
                if truth != M[idx[u]][idx[v]]:
                    bad += 1
                    if bad == 1:
                        print("   MISMATCH", p, u, v, truth, M[idx[u]][idx[v]])
    print(f"   Monk's rule vs divided-difference Schubert instrument, n={n}: "
          f"{tested-bad}/{tested} agree, {bad} mismatches")
    return bad

if __name__ == "__main__":
    for n in (3, 4, 5):
        check_monk(n)
