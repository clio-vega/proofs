"""Q212.  Of the (w,r,t) triples whose prefix set A_t has NO maximum in the
right weak order (the 178 of 2026-09-21-c1 section 9), how many have w 321-avoiding?

A_t = { w^1 w^2 ... w^t : w = w^1...w^r cyclically decreasing, length-additive }.
u <=_R v  iff  l(u) + l(u^{-1} v) = l(v).
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-q215')
sys.path.insert(0, '/home/clio/projects/proofs/code-q209')
from cylact import is_321_avoiding_words, is_321_avoiding_pattern
from affstan import identity, rmul_s, length, cyc_dec_elements, elements_of_length

def inv(w, n):
    """inverse in window notation."""
    out = [0] * n
    for i in range(1, n + 1):
        q, r = divmod(w[i - 1] - 1, n)
        out[r] = i - q * n
    return tuple(out)

def mul(a, b, n):
    """(a*b)(i) = a(b(i))."""
    def av(i):
        q, r = divmod(i - 1, n)
        return a[r] + q * n
    return tuple(av(b[i]) for i in range(n))

def leq_R(u, v, n):
    return length(u, n) + length(mul(inv(u, n), v, n), n) == length(v, n)

def prefix_sets(w, n, r):
    """Returns {t: set of achievable prefixes} for t = 1..r-1."""
    cd = cyc_dec_elements(n)
    lw = length(w, n)
    layer = {identity(n): [[identity(n)]]}       # element -> list of prefix chains
    chains = [(identity(n),)]
    for step in range(r):
        new = []
        for ch in chains:
            u = ch[-1]
            lu = length(u, n)
            for Sf, (uS, size, word) in cd.items():
                if lu + size > lw:
                    continue
                v = u
                ok = True
                for i in word:
                    v2 = rmul_s(v, i, n)
                    if length(v2, n) != length(v, n) + 1:
                        ok = False; break
                    v = v2
                if ok:
                    new.append(ch + (v,))
        chains = new
    chains = [c for c in chains if c[-1] == w]
    return {t: {c[t] for c in chains} for t in range(1, r)}

def run():
    tot = nomax = 0
    nomax_avoid = 0
    tot_avoid = 0
    examples = []
    for n in (3, 4):
        byl = elements_of_length(n, 5)
        for L in range(1, 6):
            for w in byl[L]:
                av = is_321_avoiding_words(w, n)
                for r in (2, 3, 4):
                    ps = prefix_sets(w, n, r)
                    for t, A in ps.items():
                        if len(A) <= 1:
                            continue
                        tot += 1
                        if av: tot_avoid += 1
                        has_max = any(all(leq_R(b, a, n) for b in A) for a in A)
                        if not has_max:
                            nomax += 1
                            if av:
                                nomax_avoid += 1
                                if len(examples) < 5:
                                    examples.append((n, w, r, t, sorted(A)))
    return tot, tot_avoid, nomax, nomax_avoid, examples

if __name__ == '__main__':
    tot, tot_avoid, nomax, nomax_avoid, ex = run()
    print(f"triples (w,r,t) with |A_t|>1        : {tot}   (321-avoiding w: {tot_avoid})")
    print(f"of these, A_t has NO right-weak max : {nomax}  (321-avoiding w: {nomax_avoid})")
    for e in ex:
        print("   example:", e)
