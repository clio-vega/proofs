"""
Q460 -- second engine: mod-p evaluation (independent of the sympy path).

Computes, on the tie-free stratum S=(1^{2N}):
  (i)  rank(M) over F_p  -- a certified lower bound on any state space (Lemma).
  (ii) #rows up to PROPORTIONALITY -- the sharper bound D >= #projective rows.
  (iii) the same at q^3 = -1 (a primitive 6th root of unity mod p), where the
        paper reports a determinant formula IS known: a DROP there is the
        control that rank is measuring the right thing.

Also verifies the explicit even/odd construction used in the written proof.
"""
import itertools, random, sys
P = 2147483647                      # 2^31-1, prime, = 1 mod 3 so sqrt(-3) exists

def inv(a): return pow(a % P, P-2, P)

# ---- tie-free scattering factors, mod p ---------------------------------
def dA(x, y, q): return (y-q*q%P*x) * (y-pow(q,3,P)*x) % P * inv((y-x)*(y-q*x)) % P
def dB(x, y, q): return (y-q*q%P*x) * (q*q%P*y-x) % P * inv((y-x)*(y-x)) % P
def dC(x, y, q): return pow(y-q*q%P*x, 2, P) * (x-q*y) % P * inv((y-x)*(y-x)%P*(y-q*x)) % P

def overlap(u, v): return not (u[1] < v[0] or v[1] < u[0])

def Delta(u, v, x, y, q):
    """Tie-free Delta_{U,V}(x,y) mod p.  u,v arcs with disjoint column sets."""
    u1, u2 = u; v1, v2 = v
    if u2 < v1:                 return dA(x, y, q)
    if v2 < u1:                 return dA(y, x, q)          # A^T
    if u1 < v1 and v2 < u2:     return dB(x, y, q)
    if v1 < u1 and u2 < v2:     return dB(y, x, q)          # B^T
    if u1 < v1 < u2 < v2:       return dC(x, y, q)
    if v1 < u1 < v2 < u2:       return dC(y, x, q)          # C^T
    raise ValueError("not tie-free: %r %r" % (u, v))

def matchings(cols):
    cols = sorted(cols)
    if not cols:
        yield (); return
    first, rest = cols[0], cols[1:]
    for k in range(len(rest)):
        for m in matchings(rest[:k] + rest[k+1:]):
            yield ((first, rest[k]),) + m

def ordered_matchings(cols):
    return [p for m in matchings(cols) for p in itertools.permutations(m)]

def rank_modp(rows):
    rows = [r[:] for r in rows]; n = len(rows); m = len(rows[0]) if n else 0
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if rows[i][c] % P), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        iv = inv(rows[r][c])
        rows[r] = [v*iv % P for v in rows[r]]
        for i in range(n):
            if i != r and rows[i][c] % P:
                f = rows[i][c]
                rows[i] = [(a - f*b) % P for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == n: break
    return r

def projective_rows(rows):
    """Number of rows up to scalar proportionality (normalise by 1st nonzero)."""
    sigs = set()
    for r in rows:
        k = next((i for i, v in enumerate(r) if v % P), None)
        if k is None: sigs.add(None); continue
        iv = inv(r[k])
        sigs.add(tuple(v*iv % P for v in r))
    return len(sigs)

def build(Cpre, Csuf, xs, q):
    Ps, Qs = ordered_matchings(Cpre), ordered_matchings(Csuf)
    a, b = len(Cpre)//2, len(Csuf)//2
    M = []
    for Pm in Ps:
        row = []
        for Qm in Qs:
            v = 1
            for i in range(a):
                for j in range(b):
                    v = v * Delta(Pm[i], Qm[j], xs[i], xs[a+j], q) % P
            row.append(v)
        M.append(row)
    return M, len(Ps), len(Qs)

def run(a, b, q, label, seed=7, split="evenodd"):
    random.seed(seed)
    if split == "evenodd":
        # the construction used in the written proof: prefix on EVEN columns,
        # suffix on ODD columns, fully interleaved.
        n = 2*(a+b)
        Cpre = [2*i for i in range(1, a+1)] + [2*i for i in range(a+1, a+1)]
        Cpre = [2*i for i in range(1, 2*a+1)]          # even cols 2,4,...,4a
        Csuf = [2*i-1 for i in range(1, 2*b+1)]        # odd  cols 1,3,...,4b-1
        # shift odd columns to interleave across the whole range
        Csuf = sorted(set(range(1, 4*max(a,b)+2, 2)))[:2*b]
    xs = []
    while len(xs) < a+b:
        v = random.randrange(2, P)
        if all(v != w and v != q*w % P and w != q*v % P for w in xs): xs.append(v)
    M, nP, nQ = build(Cpre, Csuf, xs, q)
    r = rank_modp(M); pr = projective_rows(M)
    print("  a=%d b=%d  %-22s M %4dx%-4d  rank=%-4d  proj-rows=%-4d  (full=%d)  %s"
          % (a, b, label, nP, nQ, r, pr, nP, "FULL" if pr == nP else "DEGENERATE"))
    return r, pr, nP

if __name__ == "__main__":
    # q with q^3 = -1: root of q^2 - q + 1 = 0, i.e. q = (1 +- sqrt(-3))/2
    s = pow(P-3, (P+1)//4, P)  # try; verify below
    q6 = None
    for cand in range(2, 100000):
        if (cand*cand - cand + 1) % P == 0: q6 = cand; break
    if q6 is None:
        # solve properly: sqrt(-3) via Tonelli for p = 3 mod 4 -> p=2^31-1 is 3 mod 4
        rt = pow((P-3) % P, (P+1)//4, P)
        assert rt*rt % P == (P-3) % P, "sqrt(-3) failed"
        q6 = (1 + rt) * inv(2) % P
    assert pow(q6, 3, P) == P-1, "q^3 != -1"
    assert (q6*q6 - q6 + 1) % P == 0
    print("q6 = %d  has q6^3 = -1 mod p  (primitive 6th root of unity)\n" % q6)

    qg = 123456789   # generic
    print("GENERIC q = %d:" % qg)
    for (a, b) in [(1,1),(2,2),(3,3),(2,3),(3,2)]:
        run(a, b, qg, "generic q")
    print("\nROOT OF UNITY q^3 = -1  (where the paper reports a determinant exists):")
    for (a, b) in [(1,1),(2,2),(3,3),(2,3),(3,2)]:
        run(a, b, q6, "q^3=-1")
