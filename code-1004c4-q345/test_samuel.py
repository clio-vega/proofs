import sys
from itertools import permutations
from collections import Counter
from schub import length, refl_length, structure_constants, rmul_t
from chains import T_tensor, T_free

def run(n, lenfn, lenname, perturb=None, verbose=False):
    C, tab = structure_constants(n, n+1)
    perms = list(permutations(range(1, n+1)))
    e = tuple(range(1, n+1))
    byl = {}
    for v in perms:
        byl.setdefault(lenfn(v), []).append(v)
    Tv = {v: T_tensor(e, v, n, lenfn) for v in perms}
    ok = bad = 0
    firstbad = None
    for u in perms:
        for w in perms:
            m = lenfn(w)-lenfn(u)
            if m < 0:
                continue
            lhs = T_tensor(u, w, n, lenfn)
            rhs = Counter()
            for v in byl.get(m, []):
                c = C.get((u, v, w), 0)
                if perturb is not None and (u, v, w) == perturb:
                    c = c+1                     # planted off-by-one (H3)
                if c:
                    for t, k in Tv[v].items():
                        rhs[t] += c*k
            lhs = Counter({k: v for k, v in lhs.items() if v})
            rhs = Counter({k: v for k, v in rhs.items() if v})
            if lhs == rhs:
                ok += 1
            else:
                bad += 1
                if firstbad is None:
                    firstbad = (u, w, lhs, rhs)
    tag = f"n={n}  length={lenname}" + ("  [PERTURBED]" if perturb else "")
    print(f"{tag:46s}  agree {ok:6d}   disagree {bad:6d}   triples(u,w) tested {ok+bad}")
    if firstbad and verbose:
        u, w, lhs, rhs = firstbad
        print("   first disagreement:  u =", u, " w =", w)
        print("     LHS T_{w/u} =", dict(lhs))
        print("     RHS sum_v c T_v =", dict(rhs))
    return ok, bad, firstbad

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    print("=== READING 1: 'strong order' = Bruhat order, length = Coxeter length ===")
    run(n, length, "Coxeter (Bruhat)")
    print()
    print("=== READING 2: 'strong order' = absolute order, length = reflection length ===")
    run(n, refl_length, "reflection (absolute)", verbose=True)
    print()
    print("=== NEGATIVE CONTROL (H3): plant +1 on one structure constant ===")
    perms = list(permutations(range(1, n+1)))
    C, _ = structure_constants(n, n+1)
    probe = sorted(C.keys())[len(C)//2]
    print("   perturbing c_{u,v}^w at", probe)
    run(n, length, "Coxeter (Bruhat)", perturb=probe, verbose=True)
