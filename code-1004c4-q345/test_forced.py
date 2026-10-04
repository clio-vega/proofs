"""(T1') chains <-> Monk matrices at n=5, small m.
   (T4)  the relation is FORCED: in degree 1 the obstruction span is exactly
         ker(W -> V), i.e. the relations (a,b)+(b,c)=(a,c)."""
from itertools import permutations, product
from fractions import Fraction
from schub import length, structure_constants
from chains import T_tensor
from monk import monk_matrix

def T1prime(n, mmax):
    perms = list(permutations(range(1, n+1)))
    idx = {w: i for i, w in enumerate(perms)}
    Ms = {p: monk_matrix(n, p)[0] for p in range(1, n)}
    bad = tested = pairs = 0
    for u in perms:
        for w in perms:
            m = length(w)-length(u)
            if m < 0 or m > mmax:
                continue
            pairs += 1
            T = T_tensor(u, w, n)
            for pp in product(range(1, n), repeat=m):
                vec = [0]*len(perms); vec[idx[u]] = 1
                for p in pp:
                    nv = [0]*len(perms)
                    for i, vi in enumerate(vec):
                        if vi:
                            for j, mm in enumerate(Ms[p][i]):
                                if mm:
                                    nv[j] += vi*mm
                    vec = nv
                tested += 1
                if vec[idx[w]] != T.get(pp, 0):
                    bad += 1
    print(f"(T1') n={n}, m<={mmax}: {pairs} (u,w) pairs, {tested} label words, "
          f"{tested-bad} agree, {bad} mismatches")

def T4(n):
    """Every (a,b) with a<b occurs as a Bruhat cover label u <| u t_{ab};
    hence the degree-1 obstructions span exactly ker(W -> V)."""
    from schub import rmul_t
    perms = list(permutations(range(1, n+1)))
    seen = set()
    for u in perms:
        lu = length(u)
        for a in range(1, n+1):
            for b in range(a+1, n+1):
                if length(rmul_t(u, a, b)) == lu+1:
                    seen.add((a, b))
    allpairs = {(a, b) for a in range(1, n+1) for b in range(a+1, n+1)}
    print(f"(T4) n={n}: reflections occurring as a cover label: "
          f"{len(seen)}/{len(allpairs)}   missing: {sorted(allpairs-seen)}")
    # and check the degree-1 obstruction for each cover equals t_ab - sum_{p=a}^{b-1} t_{p,p+1}
    C, _ = structure_constants(n, n+1)
    sp = {p: tuple(list(range(1, p))+[p+1, p]+list(range(p+2, n+1))) for p in range(1, n)}
    bad = 0; tested = 0
    for u in perms:
        lu = length(u)
        for a in range(1, n+1):
            for b in range(a+1, n+1):
                w = rmul_t(u, a, b)
                if length(w) != lu+1:
                    continue
                tested += 1
                # RHS of Samuel in degree 1, in the FREE space: sum_p c_{u,s_p}^w * t_{p,p+1}
                coeffs = {p: C.get((u, sp[p], w), 0) for p in range(1, n)}
                want = {p: (1 if a <= p < b else 0) for p in range(1, n)}
                if coeffs != want:
                    bad += 1
    print(f"(T4) n={n}: degree-1 obstruction = t_ab - sum_{{p=a}}^{{b-1}} t_{{p,p+1}} "
          f"for {tested-bad}/{tested} covers")

if __name__ == "__main__":
    T1prime(5, 4)
    for n in (3, 4, 5):
        T4(n)
