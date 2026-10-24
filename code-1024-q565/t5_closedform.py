"""T5.  (i) r=1 closed form via the kernel H(z) = prod_t (z w_t - 1)/(z - w_t).

    A_m(w) := sum_{j=0}^{min(m,n)} (-1)^j e_{n-j}(w) h_{m-j}(w)      [ = [z^-m]_inf H ]
    Psi^C_{(1),(m)}(M) = A_m / e                 (m >= 1)
    Psi^C_{(1),(0)}(M) = (e^2-1)/e^2             (m = 0)
    Psi^C_{(1),(m)}(M) = -Atil_{-m} / e          (m <= -1),  Atil = expansion at z=0

(ii) determinant ansatz for r=2 in terms of a_k := Psi^C_{(1),(k)}(M).
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
sys.path.insert(0, '/home/clio/projects/proofs/code-1024-q565')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt

def e_elem(w, j):
    n = len(w)
    if j < 0 or j > n: return F(0)
    return sum((lambda S: __import__('functools').reduce(lambda a,b: a*b, S, F(1)))(S)
               for S in itertools.combinations(w, j))

def h_comp(w, k):
    """complete homogeneous h_k(w)."""
    if k < 0: return F(0)
    if k == 0: return F(1)
    n = len(w)
    # recursion h_k = sum_{i=1..n} ... use generating: h_k = sum over multisets
    tot = F(0)
    for comb in itertools.combinations_with_replacement(range(n), k):
        p = F(1)
        for i in comb: p *= w[i]
        tot += p
    return tot

def A_m(w, m):
    n = len(w)
    return sum((-1)**j * e_elem(w, n-j) * h_comp(w, m-j) for j in range(0, min(m, n)+1))

def A_til(w, k):
    """[z^k]_0 of H(z) = (1/e) prod (1 - z w_t)/(1 - z/w_t).
       = (1/e) * sum_j (-1)^j e_j(w) h_{k-j}(1/w)."""
    e = F(1)
    for x in w: e *= x
    winv = [1/x for x in w]
    n = len(w)
    return (sum((-1)**j * e_elem(w, j) * h_comp(winv, k-j) for j in range(0, min(k, n)+1)))/e

random.seed(777)
VALS = [F(2), F(3), F(5), F(7), F(-2), F(-3), F(4), F(-5), F(3,2), F(-5,2), F(7,3)]

print("(i) r=1 closed form")
ok = bad = 0; fails = []
for n in range(1, 6):
    for trial in range(12):
        w = random.sample(VALS, n)
        if len(set(w)) != n: continue
        e = F(1)
        for x in w: e *= x
        for m in range(-3, 6):
            try: lhs = psi_headsum(w, [m])
            except ZeroDivisionError: continue
            if m >= 1:   rhs = A_m(w, m)/e
            elif m == 0: rhs = (e*e - 1)/(e*e)
            else:        rhs = -A_til(w, -m)/e
            if lhs == rhs: ok += 1
            else:
                bad += 1
                if len(fails) < 5: fails.append((n, m, w, lhs, rhs))
print(f"    Psi^C_(1),(m) = kernel coefficient / e :  {ok}/{ok+bad}   (n<=5, -3<=m<=5)")
for f in fails: print("      FAIL", f)

print()
print("(ii) r=2: is Psi^C_{(1^2),theta} a 2x2 determinant in the a_k = Psi^C_{(1),(k)} ?")
def a(w, k): return psi_headsum(w, [k])
hits = {}
for n in range(2, 6):
    for trial in range(3):
        w = random.sample(VALS, n)
        if len(set(w)) != n: continue
        for th in [(0,0),(0,1),(1,0),(1,1),(0,2),(2,0),(1,2),(2,1)]:
            m = [th[0], th[1]+1]
            if m[0] == m[1]: continue
            try: lhs = psi_headsum(w, list(th))
            except ZeroDivisionError: continue
            m1, m2 = m
            cands = {
              "a_m1*a_m2":                  a(w,m1)*a(w,m2),
              "det[[a_m1,a_m2],[a_m1+1,a_m2+1]]": a(w,m1)*a(w,m2+1)-a(w,m2)*a(w,m1+1),
              "det[[a_m1,a_m1+1],[a_m2,a_m2+1]]": a(w,m1)*a(w,m2+1)-a(w,m1+1)*a(w,m2),
              "a_m1*a_m2 - a_m1+1*a_m2-1":  a(w,m1)*a(w,m2)-a(w,m1+1)*a(w,m2-1),
            }
            for name, v in cands.items():
                hits.setdefault(name, [0,0])
                if v == lhs: hits[name][0] += 1
                else: hits[name][1] += 1
for name, (g, b) in hits.items():
    print(f"    {name:38} {g}/{g+b}")
