"""The layer profile at general m.

Each f_nu is a convolution of m interval indicators, symmetric about the common
centre c=(d-b)/2 (concentricity).  Write gamma_nu(y)=f_nu(c+y) and
beta_nu(y)=gamma_nu(y)-gamma_nu(y+1) >= 0: the LAYER PROFILE, i.e.
f_nu = sum_r beta_nu(r) * 1_[c-r, c+r].
Then G = sum_nu f_nu has gamma(y)=sum_r>=y beta(r) with beta = sum_nu beta_nu.
By the tail lemma (2026-10-01, proved):  beta log-concave  ==>  G log-concave.

At m=2, beta was the POSITIVE PART OF A CONCAVE function.  Question: at m>=3?
"""
import gen
from itertools import product

def spline(widths):
    return gen.conv_intervals(list(widths))

def radial(co):
    """co: symmetric coeff list.  Returns gamma as list indexed by half-steps from centre:
    gamma[j] = co[centre_index_outward j].  Works for odd and even length."""
    L = len(co)
    if L == 0: return []
    # centre index (L-1)/2 ; distances are 0,1,2,... if L odd, 1/2,3/2,... if L even
    if L % 2 == 1:
        mid = (L - 1) // 2
        return [co[mid + j] for j in range(mid + 1)]
    else:
        mid = L // 2
        return [co[mid + j] for j in range(mid)]

def beta_of(gam):
    return [gam[j] - (gam[j+1] if j+1 < len(gam) else 0) for j in range(len(gam))]

def is_lc(v):
    """log-concave with interval support, on a list indexed 0.. (index 0 = innermost layer)"""
    for i in range(len(v)):
        l = v[i-1] if i-1 >= 0 else 0
        r = v[i+1] if i+1 < len(v) else 0
        if v[i]*v[i] < l*r: return False
    nz=[i for i,c in enumerate(v) if c]
    if nz and nz[-1]-nz[0]+1 != len(nz): return False
    return True

def is_concave_pospart(v):
    """is v = c_+ for some concave integer c?  Equivalent, for v>=0 with interval
    support [p,q]: v is concave on [p-1, q+1] when extended by 0 outside, i.e.
    v concave on its support AND the boundary slopes admit extension.
    Simplest exact test: v is the positive part of a concave sequence iff
    v restricted to its support is concave and its first/last increments allow
    a concave continuation to 0, which is automatic.  So: concave on support."""
    nz=[i for i,c in enumerate(v) if c]
    if not nz: return True
    if nz[-1]-nz[0]+1 != len(nz): return False
    s=v[nz[0]:nz[-1]+1]
    inc=[s[i+1]-s[i] for i in range(len(s)-1)]
    return all(inc[i] >= inc[i+1] for i in range(len(inc)-1))

if __name__ == "__main__":
    # ---- Q1: is beta_w log-concave for a SINGLE width vector w? ----
    print("Q1: single-spline layer profile beta_w")
    bad_lc=[]; bad_cc=[]; tot=0
    for m in range(1, 6):
        for w in product(range(1, 9), repeat=m):
            tot+=1
            b=beta_of(radial(spline(w)))
            if not is_lc(b): bad_lc.append((w,b))
            if not is_concave_pospart(b): bad_cc.append((w,b))
    print(f"  tested {tot} width vectors, m<=5, w_i<=8")
    print(f"  beta_w NOT log-concave: {len(bad_lc)}   first: {bad_lc[:3]}")
    print(f"  beta_w NOT positive part of concave: {len(bad_cc)}   first: {bad_cc[:4]}")
