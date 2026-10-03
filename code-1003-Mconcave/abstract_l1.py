"""THE ABSTRACT CLASS TEST (brief hazard 3).

Target in reduced form.  B = prod_i [P_i,Q_i] a box in Z^m, sigma in Z.
On the slice {sum y_i = sigma} we have ||y||_1 = sigma + 2*k(y) with
k(y) = sum_i max(-y_i,0).  Since Lambda = (n-m)/2 - ||nu-a||_1/2, the half-width
profile M(r) is (up to r = const - k) the function

    N(k) = #{ y in B : sum y_i = sigma, sum_i max(-y_i,0) = k }.

QUESTION: is N always the positive part of a concave function, for an ARBITRARY
box and arbitrary sigma?  If yes -> the theorem is a pure lattice-point statement
and real slices are irrelevant.  If no -> I need the real constraints on P_i,Q_i.
"""
from itertools import product
from collections import Counter
from layers import is_concave_pospart

def profile(P, Q, sigma):
    m = len(P)
    c = Counter()
    for y in product(*[range(P[i], Q[i]+1) for i in range(m)]):
        if sum(y) != sigma: continue
        c[sum(max(-t,0) for t in y)] += 1
    if not c: return []
    return [c.get(k,0) for k in range(min(c), max(c)+1)]

st = Counter(); wit = []
RNG = range(-3, 4)
for m in (2,3,4):
    for P in product(RNG, repeat=m):
        for Q in product(RNG, repeat=m):
            if any(Q[i] < P[i] for i in range(m)): continue
            for sigma in range(sum(P), sum(Q)+1):
                N = profile(P,Q,sigma)
                if not N: continue
                st[f'm{m}_tot'] += 1
                if is_concave_pospart(N): st[f'm{m}_ok'] += 1
                else:
                    st[f'm{m}_FAIL'] += 1
                    if len(wit) < 12: wit.append((m,P,Q,sigma,N))
print(dict(st))
print("failures:")
for w in wit: print('  ', w)
