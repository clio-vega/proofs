"""Route observables for Q215.  Nothing here tests M-convexity of F~_w (that is
a theorem of WZZ and therefore unfalsifiable -- see PROVE.md section 3).  What is
tested is the DICTIONARY and the two new lemmas:

  T0  u_i is a representation: u_w . mu is independent of the reduced word.
  T1  For 321-avoiding w there is a cylindric shape mu with u_w . mu != 0
      (Lam thm:321), and for NON-321-avoiding w there is none (Lam's converse).
  T2  (D1) COEFFICIENTWISE:  F~_w(x_1..x_r) = s^c_{lam/mu}(x_1..x_r) for all r,
      where both sides are computed from their own definitions.  This is the
      dictionary test and it can fail.
  T3  (D2) lambda-hat(lam/mu, ell) is INDEPENDENT of ell for ell >= ellmin.
  T4  (Q216) lambda-hat = mu(w), Lam thm:monomial's conjugate-of-sorted-code.
      Two unrelated routes to the same partition; untuned.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-q215')
sys.path.insert(0, '/home/clio/projects/proofs/code-q209')
from cylact import *
from affstan import elements_of_length, affine_stanley_support, length

def all_shapes(n, m):
    """Canonical window representatives: 0 <= x_1 < ... < x_m < n."""
    return [tuple(c) for c in itertools.combinations(range(n), m)]

def realise(w, n, maxlen=10):
    """All (m, mu, lam) with u_w . mu = lam != 0.  Also checks reduced-word
    independence.  Returns (list, n_words)."""
    words = reduced_words(w, n, maxlen=maxlen)
    out = []
    for m in range(1, n):
        for mu in all_shapes(n, m):
            vals = {u_word(mu, wd, n, m) for wd in words}
            if len(vals) != 1:
                raise AssertionError(f"reduced-word dependence: {w} {mu}")
            lam = vals.pop()
            if lam is not None:
                out.append((m, mu, lam))
    return out, len(words)

LR = {3: (9, 6), 4: (7, 5), 5: (6, 4), 6: (5, 4), 7: (4, 3)}

def run(nmax=5, report=print):
    stats = dict(T0=0, T1a=0, T1b=0, T1b_fail=0, T2=0, T2_fail=0,
                 T3=0, T3_fail=0, T4=0, T4_fail=0, realised=0, unrealised=0)
    fails = []
    for n in range(3, nmax + 1):
        L, R = LR[n]
        byl = elements_of_length(n, L)
        for l in range(1, L + 1):
            for w in byl[l]:
                av = is_321_avoiding_words(w, n)
                reals, nw = realise(w, n)
                stats['T0'] += nw
                if av:
                    stats['T1a'] += 1
                    if not reals:
                        fails.append(('T1a no realisation', n, w)); stats['unrealised'] += 1
                        continue
                    stats['realised'] += 1
                else:
                    stats['T1b'] += 1
                    if reals:
                        stats['T1b_fail'] += 1
                        fails.append(('T1b non-avoiding realised', n, w, reals[:2]))
                    continue
                # --- T4: lambda-hat vs mu(w)
                muw = mu_of_w(w, n)
                for (m, mu, lam) in reals:
                    lh, ellmin = lambda_hat(mu, lam, n, m)
                    stats['T4'] += 1
                    if lh != muw:
                        stats['T4_fail'] += 1
                        fails.append(('T4', n, w, m, mu, lam, lh, muw))
                    # --- T3: lambda-hat independent of ell
                    gammas, em = greedy_chain(mu, lam, n, m)
                    for ell in range(em, em + 4):
                        stats['T3'] += 1
                        lh2 = tuple(sorted((x for x in (gammas + [0]*4)[:ell] if x > 0),
                                           reverse=True))
                        if lh2 != lh:
                            stats['T3_fail'] += 1
                            fails.append(('T3', n, w, m, mu, lam, ell, lh2, lh))
                    # --- T2: coefficientwise identity in r variables
                    for r in range(1, R + 1):
                        A = affine_stanley_support(w, n, r)
                        B = cyl_poly_support(mu, lam, n, m, r)
                        stats['T2'] += 1
                        if A != B:
                            stats['T2_fail'] += 1
                            if len(fails) < 40:
                                fails.append(('T2', n, w, m, mu, lam, r,
                                              sorted(A.items())[:6], sorted(B.items())[:6]))
    return stats, fails

if __name__ == '__main__':
    s, f = run(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
    for k, v in s.items():
        print(f"{k:12s} {v}")
    print(f"total failures: {len(f)}")
    for x in f[:15]:
        print("  ", x)
