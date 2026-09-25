"""Q254 calibration: is my transcription of FMS's chi_D correct?

Independent check.  FMS Theorem 4 (l.239, cited to Kraskiewicz-Pragacz) says
chi_{D(w)} = Schubert polynomial S_w; FMS Theorem 5 (l.243, cited to
\bibitem{keypolynomials} = Reiner-Shimozono, JCTA 70 (1995) 107-143 -- NB the
.bbl has the titles of \bibitem{demazure} and \bibitem{keypolynomials} swapped)
says chi_{D(alpha)} = key polynomial kappa_alpha.

Here Schubert and key polynomials are computed FROM SCRATCH by divided
differences -- no matroids, no flagged Weyl modules -- and their supports are
compared with supp(chi_D) computed by fms.supp_chi.  A disagreement would mean
I mis-transcribed a definition.
"""
import sys, itertools
from collections import defaultdict
sys.path.insert(0, '/home/clio/projects/proofs/code-q254')
from fms import supp_chi

# ---- polynomials as dict: exponent tuple -> coefficient ----------------
def padd(f, g):
    h = defaultdict(int)
    for k, v in f.items(): h[k] += v
    for k, v in g.items(): h[k] += v
    return {k: v for k, v in h.items() if v}

def pmulvar(f, i, n):
    out = {}
    for k, v in f.items():
        e = list(k); e[i] += 1
        out[tuple(e)] = v
    return out

def swap(f, i):
    out = {}
    for k, v in f.items():
        e = list(k); e[i], e[i+1] = e[i+1], e[i]
        out[tuple(e)] = out.get(tuple(e), 0) + v
    return out

def divdiff(f, i, n):
    """(f - s_i f)/(x_i - x_{i+1}).  Done monomial-pair-wise:
       (x_i^a x_{i+1}^b - x_i^b x_{i+1}^a)/(x_i-x_{i+1}) = sum_{c=b}^{a-1} x_i^c x_{i+1}^{a+b-1-c}
       for a>b, = 0 for a=b, = -(...) for a<b."""
    num = padd(f, {k: -v for k, v in swap(f, i).items()})
    out = defaultdict(int)
    done = set()
    for k, v in num.items():
        if k in done: continue
        a, b = k[i], k[i+1]
        if a <= b: continue
        kk = list(k); kk[i], kk[i+1] = b, a; kk = tuple(kk)
        done.add(k); done.add(kk)
        # coefficient of x^k in num is v; of x^kk is num.get(kk)
        for c in range(b, a):
            e = list(k); e[i] = c; e[i+1] = a + b - 1 - c
            out[tuple(e)] += v
    return {k: v for k, v in out.items() if v}

def schubert(w, n):
    """Schubert polynomial of w in S_n, by FMS l.78-88's recursion."""
    w = tuple(w)
    w0 = tuple(range(n, 0, -1))
    if w == w0:
        return {tuple(n - 1 - i for i in range(n)): 1}
    for i in range(n - 1):
        if w[i] < w[i+1]:
            ws = list(w); ws[i], ws[i+1] = ws[i+1], ws[i]
            return divdiff(schubert(tuple(ws), n), i, n)
    raise RuntimeError

def key(alpha, n):
    """Key polynomial kappa_alpha, FMS l.92-96."""
    alpha = tuple(alpha)
    if all(alpha[i] >= alpha[i+1] for i in range(len(alpha)-1)):
        return {alpha: 1}
    for i in range(len(alpha)-1):
        if alpha[i] < alpha[i+1]:
            ah = list(alpha); ah[i], ah[i+1] = ah[i+1], ah[i]
            return divdiff(pmulvar(key(tuple(ah), n), i, n), i, n)
    raise RuntimeError

def rothe(w, n):
    """D(w)_j = {i : w(i)>j, w^{-1}(j)>i},  FMS l.97-99."""
    winv = [0]*(n+1)
    for i, wi in enumerate(w, 1): winv[wi] = i
    return [tuple(i for i in range(1, n+1) if w[i-1] > j and winv[j] > i)
            for j in range(1, n+1)]

def skyline(alpha, n):
    """D(alpha)_j = {i : alpha_i >= j}.  (FMS l.147 PRINTS '{j<=n: alpha_j>=j}',
    which is a typo -- their own Figure 2, alpha=(3,2,1,0,1) -> ({1,2,3,5},...),
    requires {i : alpha_i >= j}.)"""
    return [tuple(i for i in range(1, n+1) if (alpha[i-1] if i-1 < len(alpha) else 0) >= j)
            for j in range(1, n+1)]

if __name__ == "__main__":
    print("=== FMS Thm 4 (Rothe diagram -> Schubert polynomial): supports ===")
    for n in range(2, 6):
        bad = 0; tot = 0
        for w in itertools.permutations(range(1, n+1)):
            tot += 1
            S1 = set(schubert(w, n).keys())
            S2 = supp_chi(rothe(w, n), n)
            if S1 != S2:
                bad += 1
                if bad <= 3: print("   MISMATCH", w, sorted(S1), sorted(S2))
        print(f"   n={n}: {tot-bad}/{tot} agree")

    print("=== FMS Thm 5 (skyline diagram -> key polynomial): supports ===")
    for n in range(2, 5):
        bad = 0; tot = 0
        for alpha in itertools.product(range(0, n+1), repeat=n):
            tot += 1
            S1 = set(key(alpha, n).keys())
            S2 = supp_chi(skyline(alpha, n), n)
            if S1 != S2:
                bad += 1
                if bad <= 3: print("   MISMATCH", alpha, sorted(S1), sorted(S2))
        print(f"   n={n}: {tot-bad}/{tot} agree")
