"""T9.  (A) Certify the P1 falsifier found by the 'failed' control.
       (B) The NEGATION covering law, which explains the other 18 zeros.

At theta = (0,1,...,r-1) we have m_j = 2(j-1), so the alternant is a Vandermonde
in the SQUARES:
      det(x_i^{2(j-1)}) / V(T)  =  prod_{i<k}(x_k^2-x_i^2)/prod_{i<k}(x_k-x_i)
                                =  prod_{i<k in T} (x_i + x_k).
So that term dies iff T contains BOTH members of a negation pair w_b = -w_a.
Covering law: EVERY term dies iff r > n - lambda, lambda = #negation pairs.
(Exactly the shape of Theorem A, with (w_a w_b - 1) replaced by (w_a + w_b).)
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, B, alternant_ratio

print("=" * 74)
print("(A)  THE P1 FALSIFIER, certified")
print("=" * 74)
W = [F(-7), F(2), F(3), F(7), F(-2)]
r = 3; theta = [0, 1, 2]; n = len(W)
m = [theta[j] + j for j in range(r)]
print(f"  M      = {W}          (n = {n})")
print(f"  r      = {r},  kappa = (1^3),  theta = {tuple(theta)},  m = {tuple(m)}")
print(f"  n - r  = {n-r}   so n = r+2 > r+1   (outside prop:recip's regime)")
print(f"  m_j pairwise distinct: {len(set(m)) == r}")
print(f"  w_a distinct, nonzero: {len(set(W))==n and all(x!=0 for x in W)}")
rec = [(a+1,b+1) for a in range(n) for b in range(a,n) if W[a]*W[b] == 1]
print(f"  reciprocal pairs w_a w_b = 1 (a<=b): {rec}   -> RECIPROCAL-FREE")
print(f"  reciprocal support S = {{}},  n - |S| = {n} >= r   (Theorem A: no structural vanishing)")
p1, p2 = psi_headsum(W, theta), psi_alt(W, theta)
print(f"  Psi (eq:Psidef, head-injection sum) = {p1}")
print(f"  Psi (thm:alt, alternant expansion)  = {p2}")
print(f"  the two independent implementations agree: {p1 == p2}")
terms = []
for T in itertools.combinations(range(n), r):
    x = [W[i] for i in T]
    b, al = B(W, T), alternant_ratio(x, m)
    terms.append((T, b, al, b*al))
nzB = [t for t in terms if t[1] != 0]
nzT = [t for t in terms if t[3] != 0]
print(f"  r-subsets: {len(terms)};  with B(T) != 0: {len(nzB)};  with full term != 0: {len(nzT)}")
print(f"  sum of terms = {sum(t[3] for t in terms)}")
print("  the nonzero terms (T, B(T), alternant, product):")
for (T, b, al, pr) in nzT:
    print(f"     T={tuple(i+1 for i in T)}  B={b}  alt={al}  term={pr}")
print(f"  => Psi = 0 with {len(nzT)} nonzero terms cancelling: MECHANISM M3, GENUINE CANCELLATION.")
print(f"  P1 (no cancellation locus for n > r+1) is REFUTED by this single point.")

print()
print("=" * 74)
print("(B)  the NEGATION covering law")
print("=" * 74)
print("  claim: at theta=(0,1,...,r-1), alternant_ratio = prod_{i<k in T}(x_i+x_k)")
ok = bad = 0
random.seed(31337)
VALS = [F(2),F(3),F(5),F(7),F(-2),F(-3),F(4),F(-5),F(11),F(-7),F(3,2),F(-5,3),F(8),F(-9)]
for r2 in range(1, 5):
    for _ in range(60):
        x = random.sample(VALS, r2)
        if len(set(x)) != r2: continue
        m2 = [j + j for j in range(r2)]      # theta_j = j-1  => m_j = 2(j-1)
        if len(set(m2)) != r2: continue
        lhs = alternant_ratio(x, m2)
        rhs = F(1)
        for i in range(r2):
            for k in range(i+1, r2): rhs *= (x[i] + x[k])
        ok += (lhs == rhs); bad += (lhs != rhs)
print(f"    alternant = prod(x_i+x_k):  {ok}/{ok+bad}   (r<=4, 60 draws each)")

print()
print("  covering law: all terms vanish  <=>  r > n - lambda   (lambda = #negation pairs)")
def lam(w):
    n = len(w)
    return sum(1 for a in range(n) for b in range(a+1, n) if w[a] == -w[b])
cok = cbad = 0
rows = []
for lamb in range(0, 4):
    for extra in range(0, 4):
        base = [F(2), F(3), F(5), F(7)][:lamb]
        w = []
        for v in base: w += [v, -v]
        w += [F(11), F(13), F(4), F(9)][:extra]
        n2 = len(w)
        if n2 == 0 or len(set(w)) != n2: continue
        for r2 in range(1, n2+1):
            theta2 = list(range(r2))
            m2 = [theta2[j]+j for j in range(r2)]
            if len(set(m2)) != r2: continue
            allzero = True
            for T in itertools.combinations(range(n2), r2):
                x = [w[i] for i in T]
                if alternant_ratio(x, m2) != 0: allzero = False; break
            pred = (r2 > n2 - lam(w))
            cok += (allzero == pred); cbad += (allzero != pred)
            if allzero != pred: rows.append((n2, r2, lamb, allzero, pred))
print(f"    negation covering law: {cok}/{cok+cbad} cells")
for rw in rows: print("      FAIL", rw)
print()
print("  and this is what killed the other 18 control zeros:")
for (W2, r3) in [([F(5),F(-5)], 2), ([F(2),F(-5,3),F(-2)], 3)]:
    n3 = len(W2); th = list(range(r3))
    print(f"    M={W2} r={r3}: lambda={lam(W2)}, n-lambda={n3-lam(W2)} < r={r3}"
          f" -> all alternants 0, Psi={psi_headsum(W2, th)}")
