"""T5. The reduced identity, symbolically.

   Ztil(L,s) = sum_{eps in {+-1}^L} eta^{N-1} prod_{alpha<beta, at least one free} (e_a s_a + e_b s_b)^{-1}
   N = L+s, free = 1..L (flipped), frozen = c_1..c_s (not flipped), eta = prod eps.

CLAIM:  Ztil = 0  <=>  L >= s+2.

Also verifies the two inputs of the proof:
  (a) P := Ztil * D is a POLYNOMIAL divisible by prod_{a<b<=L}(s_a^2 - s_b^2),
      where D = prod_{a<b<=L}(s_a^2-s_b^2) * prod_{a<=L,p}(c_p^2-s_a^2);
  (b) deg_{s_1} P <= N-1.
"""
import itertools, sympy as sp

def Ztil(L, s):
    sv = sp.symbols(f's1:{L+1}')
    cv = sp.symbols(f'c1:{s+1}')
    N = L+s
    vals = list(sv)+list(cv)
    tot = 0
    for eps in itertools.product([1,-1], repeat=L):
        e = list(eps)+[1]*s
        eta = sp.prod(eps)
        den = 1
        for a in range(N):
            for b in range(a+1, N):
                if a < L or b < L:          # at least one free endpoint
                    den *= (e[a]*vals[a] + e[b]*vals[b])
        tot += eta**(N-1)/den
    return sp.together(sp.simplify(sp.cancel(tot))), sv, cv

print(f"{'L':>2} {'s':>2} {'N':>2} | {'L>=s+2':>6} | Ztil")
for L in range(1, 5):
    for s in range(0, 4):
        if L+s > 6: continue
        Z, sv, cv = Ztil(L, s)
        z = sp.simplify(Z)
        isz = (z == 0)
        pred = L >= s+2
        mark = "" if isz == pred else "   <<< DISAGREES"
        txt = "0" if isz else ("nonzero  " + str(sp.factor(sp.numer(sp.cancel(z))))[:70])
        print(f"{L:>2} {s:>2} {L+s:>2} | {str(pred):>6} | {txt}{mark}")

print()
print("--- proof inputs (a) divisibility and (b) degree, on NON-vanishing cells (L<=s+1) ---")
print("    (if they held only where Ztil=0 they would be vacuous; a positive control)")
for (L, s) in [(1,1),(1,2),(2,1),(2,2),(2,3),(3,2),(3,3)]:
    if L+s > 6: continue
    Z, sv, cv = Ztil(L, s)
    N = L+s
    D = 1
    for a in range(L):
        for b in range(a+1, L): D *= (sv[a]**2 - sv[b]**2)
    for a in range(L):
        for p in range(s): D *= (cv[p]**2 - sv[a]**2)
    P = sp.cancel(sp.expand(Z*D))
    ispoly = P.is_polynomial(*(list(sv)+list(cv)))
    div = 1
    for a in range(L):
        for b in range(a+1, L): div *= (sv[a]**2 - sv[b]**2)
    q, rem = sp.div(sp.expand(P), sp.expand(div), sv[0]) if L >= 2 else (P, 0)
    divides = (sp.simplify(sp.cancel(P/div)).is_polynomial(*(list(sv)+list(cv)))) if L >= 2 else True
    deg1 = sp.Poly(sp.expand(P), sv[0]).degree() if P != 0 else -1
    print(f"  L={L} s={s} N={N}: P polynomial={ispoly}  prod(s_a^2-s_b^2) | P : {divides}"
          f"   deg_s1 P = {deg1}  (bound N-1 = {N-1}; divisor deg_s1 = {2*(L-1)})")
