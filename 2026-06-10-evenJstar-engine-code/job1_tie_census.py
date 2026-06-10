"""
JOB 1 — extend the even-|J*| tie census for G_lambda(i)=0  (general lambda, d=4).

Machinery (exact, Z[i]; pi = 1+i, v_pi(N) = 2 v_2(N) for N in Z):
    G_lambda(i) = sum_{j=0}^m C(m,j) i^j (1+i)^j M_j,
    M_j = <s_lambda, h1^{2(m-j)} e2^j>  in  Z_{>=0},   M_0 = f^lambda.
    val(j) = j + 2 v_2( C(m,j) M_j ),   J* = argmin_j val(j),   mu = min val.

Key shortcut (no Sage needed):  h1 = p1, e2 = (p1^2 - p2)/2, so
    M_j = 2^{-j} sum_{b=0}^j C(j,b)(-1)^b chi^lambda(2^b 1^{2m-2b}),
where chi^lambda(2^b 1^{2m-2b}) is a Murnaghan-Nakayama character on a
class with b parts of size 2 and the rest 1.

We compute chi_b := chi^lambda(2^b 1^{2m-2b}) for b=0..m once per shape, then
binomial-transform to all M_j.  Cross-checked against the symfunc.py p-basis
scalar product and against the documented (2,2) case M=[2,1,1].
"""
import sys, math, csv
from functools import lru_cache
from collections import Counter

# ---- partitions ----
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

# ---- Murnaghan-Nakayama (border-strip recursion) ----
@lru_cache(maxsize=None)
def mn_character(lam, mu):
    """chi^lam(mu), integer, via MN recursion. lam,mu tuples desc."""
    lam = tuple(x for x in lam if x > 0)
    mu = tuple(x for x in mu if x > 0)
    if sum(lam) != sum(mu):
        return 0
    if len(mu) == 0:
        return 1 if len(lam) == 0 else 0
    r = mu[0]
    rest = mu[1:]
    total = 0
    for (newlam, height) in border_strips(lam, r):
        total += (-1)**height * mn_character(newlam, rest)
    return total

def border_strips(lam, r):
    lam = list(lam)
    k = len(lam)
    beta = [lam[i] + (k - 1 - i) for i in range(k)]
    beta_set = set(beta)
    out = []
    for b in beta:
        if b - r >= 0 and (b - r) not in beta_set:
            new_beta = sorted([x for x in beta if x != b] + [b - r], reverse=True)
            height = sum(1 for x in beta if (b - r) < x < b)
            newlam = tuple(new_beta[i] - (k - 1 - i) for i in range(k))
            newlam = tuple(x for x in newlam if x > 0)
            out.append((newlam, height))
    return out

def v2(n):
    """2-adic valuation of a nonzero integer; +inf signalled by None for n==0."""
    if n == 0:
        return None
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c

def chi_b_vector(lam, m):
    """chi_b = chi^lam(2^b 1^{2m-2b}) for b=0..m."""
    return [mn_character(lam, tuple([2]*b + [1]*(2*m - 2*b))) for b in range(m+1)]

def M_vector(lam, m):
    """M_j for j=0..m via binomial transform of chi_b. Returns exact ints."""
    chi = chi_b_vector(lam, m)
    Ms = []
    for j in range(m+1):
        # M_j * 2^j = sum_{b=0}^j C(j,b)(-1)^b chi_b
        s = 0
        for b in range(j+1):
            s += math.comb(j, b) * ((-1)**b) * chi[b]
        assert s % (2**j) == 0, f"M_j not integer: lam={lam} j={j} s={s}"
        Ms.append(s // (2**j))
    return Ms

def G_gaussian(lam, m, Ms):
    """Exact G_lambda(i) = sum_j C(m,j) i^j (1+i)^j M_j as (re, im) integers."""
    re, im = 0, 0
    for j in range(m+1):
        coeff = math.comb(m, j) * Ms[j]
        # i^j
        ij = [(1,0),(0,1),(-1,0),(0,-1)][j % 4]
        # (1+i)^j : compute exactly
        pr, pi = 1, 0
        for _ in range(j):
            pr, pi = pr - pi, pr + pi  # *(1+i)
        # term = coeff * i^j * (1+i)^j
        ar, ai = ij
        # (ar+ai i)(pr+pi i)
        tr = ar*pr - ai*pi
        ti = ar*pi + ai*pr
        re += coeff * tr
        im += coeff * ti
    return re, im

def vpi(re, im):
    """v_pi of a Gaussian integer (re+im i), pi=1+i.  None if zero (infinite)."""
    if re == 0 and im == 0:
        return None
    c = 0
    while True:
        if (re + im) % 2 != 0:
            return c
        # divide by (1+i): (re+im i)(1-i)/2 = ((re+im)+(im-re)i)/2
        nre = (re + im) // 2
        nim = (im - re) // 2
        re, im = nre, nim
        c += 1

def analyze(lam, m):
    Ms = M_vector(lam, m)
    val = []
    for j in range(m+1):
        c = math.comb(m, j) * Ms[j]
        vv = v2(c)
        if vv is None:
            val.append(None)  # M_j = 0 -> infinite val, never the min
        else:
            val.append(j + 2*vv)
    finite = [v for v in val if v is not None]
    mu = min(finite)
    Jstar = [j for j, v in enumerate(val) if v == mu]
    re, im = G_gaussian(lam, m, Ms)
    vp = vpi(re, im)
    return dict(lam=lam, m=m, M=Ms, val=val, Jstar=Jstar, mu=mu,
               vpi=vp, re=re, im=im)

# ---------- cross-checks ----------
def selfcheck():
    # (2,2): M=[2,1,1], val=[2,3,2], J*={0,2}, G=0
    r = analyze((2,2), 2)
    assert r['M'] == [2,1,1], r['M']
    assert r['val'] == [2,3,2], r['val']
    assert r['Jstar'] == [0,2], r['Jstar']
    assert r['vpi'] is None, r  # G=0
    # cross-check M_j against p-basis scalar product (symfunc) for a few shapes
    import symfunc as sf
    for lam, m in [((4,),2), ((3,1),2), ((2,1,1),2), ((4,2),3), ((3,3),3), ((2,2,2),3)]:
        Ms = M_vector(lam, m)
        for j in range(m+1):
            e2j = sf.power(sf.e_n(2), j)
            h1pow = sf.power(sf.h1(), 2*(m-j))
            f = sf.mul(h1pow, e2j)
            ip = int(sf.inner(sf.schur_p(lam), f))
            assert ip == Ms[j], f"MISMATCH lam={lam} j={j}: shortcut={Ms[j]} pbasis={ip}"
    print("selfcheck PASSED (shortcut == p-basis; (2,2) reference matches)")

if __name__ == "__main__":
    selfcheck()
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    rows = []
    tie_rows = []
    zero_outliers = []
    jstar_size_dist = Counter()
    for m in range(1, MMAX+1):
        n = 2*m
        cnt = 0; ties = 0
        for lam in partitions(n):
            r = analyze(lam, m)
            cnt += 1
            sz = len(r['Jstar'])
            row = dict(m=m, lam="|".join(map(str,lam)), f=r['M'][0],
                       Jstar="|".join(map(str,r['Jstar'])), Jsize=sz,
                       mu=r['mu'], vpi=("inf" if r['vpi'] is None else r['vpi']),
                       parity=(r['Jstar'][0]%2))
            rows.append(row)
            jstar_size_dist[sz] += 1
            if r['vpi'] is None and lam != (2,2):
                zero_outliers.append(lam)
            if sz >= 2:
                ties += 1
                tie_rows.append(dict(m=m, lam="|".join(map(str,lam)),
                    parity=r['Jstar'][0]%2,
                    Jstar="|".join(map(str,r['Jstar'])), Jsize=sz,
                    mu=r['mu'],
                    vpi=("inf" if r['vpi'] is None else r['vpi']),
                    vpi_minus_mu=("inf" if r['vpi'] is None else r['vpi']-r['mu']),
                    Mvec="|".join(map(str,r['M'])),
                    valvec="|".join("inf" if v is None else str(v) for v in r['val'])))
        # even-|J*| conjecture: every TIE (|J*|>=2) has EVEN size.
        # (unique-min |J*|=1 is odd but handled directly by the Newton polygon.)
        bad = [r for r in rows if r['m']==m and r['Jsize']>=2 and r['Jsize']%2==1]
        uniq = sum(1 for r in rows if r['m']==m and r['Jsize']==1)
        print(f"m={m:2d}  shapes={cnt:5d}  unique-min={uniq:5d}  ties(|J*|>=2)={ties:4d}  "
              f"ODD-tie violations={len(bad)}")
    print()
    print("|J*| size distribution (all shapes, all m):", dict(sorted(jstar_size_dist.items())))
    print("G=0 outliers (lam != (2,2)):", zero_outliers if zero_outliers else "NONE")

    with open("results/tie-shapes.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["m","lam","parity","Jstar","Jsize","mu","vpi","vpi_minus_mu","Mvec","valvec"])
        w.writeheader(); w.writerows(tie_rows)
    with open("results/jstar-census.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["m","lam","f","Jstar","Jsize","mu","vpi","parity"])
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote results/tie-shapes.csv ({len(tie_rows)} ties), results/jstar-census.csv ({len(rows)} shapes)")
