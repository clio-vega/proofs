"""Residual question for Gap A.
Decomposition (PROVE.md):
  v_pi(G_lambda) = v2(f^lambda) + [ v_pi(G_4core) - v2(f^4core) ] + residual(lambda)
=> residual(lambda) = v_pi(G_lam) - v2(f^lam) - v_pi(G_4core) + v2(f^4core)
Test: is residual(lambda) a function of the 4-quotient alone?
"""
import sys, os
sys.path.insert(0, os.path.join("/home/clio/projects/scratch/2026-06-09-d4-involution"))
from core import psi_power_schur
from math import comb
from functools import reduce
from collections import defaultdict

def v2(n):
    if n == 0: return float('inf')
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

def vpi(re, im):
    # v_pi(z) = v2(N(z)) = v2(re^2+im^2); inf if zero
    if re == 0 and im == 0: return float('inf')
    return v2(re*re + im*im)

def hook_f(lam):
    """number of SYT via hook length formula (exact integer)."""
    lam = [x for x in lam if x>0]
    if not lam: return 1
    n = sum(lam)
    conj = []
    for c in range(lam[0]):
        conj.append(sum(1 for r in lam if r> c))
    prod = 1
    for i,row in enumerate(lam):
        for j in range(row):
            arm = row - (j+1)
            leg = conj[j] - (i+1)
            prod *= (arm+leg+1)
    from math import factorial
    return factorial(n)//prod

def beta_set(lam, r):
    """beta numbers with r beads (r>=len(lam)). b_i = lam_i + (r-i), i=1..r."""
    lam = list(lam) + [0]*(r-len(lam))
    return [lam[i] + (r-1-i) for i in range(r)]  # i from 0: lam_i + (r-1-i)

def beta_to_part(betas):
    betas = sorted(betas, reverse=True)
    r = len(betas)
    part = [betas[i] - (r-1-i) for i in range(r)]
    return tuple(x for x in part if x>0)

def core_quotient_4(lam):
    lam = [x for x in lam if x>0]
    ell = len(lam)
    # choose r multiple of 4, r >= ell
    r = ell
    while r % 4 != 0: r += 1
    if r == 0: r = 4
    betas = beta_set(lam, r)
    runners = {j: [] for j in range(4)}
    for b in betas:
        runners[b%4].append(b//4)  # normalized position on runner
    # quotient: each runner's normalized positions -> partition (beta-set on runner)
    quotient = []
    for j in range(4):
        positions = runners[j]  # these are the "beta numbers" on runner j
        quotient.append(beta_to_part(positions))
    # core: push beads up each runner: t beads -> positions 0..t-1
    core_betas = []
    for j in range(4):
        t = len(runners[j])
        for k in range(t):
            core_betas.append(4*k + j)
    core = beta_to_part(core_betas)
    return core, tuple(quotient)

# ---- self check ----
def selfcheck():
    import itertools
    def parts(n):
        # generate partitions of n
        def gen(n, mx):
            if n==0:
                yield ()
                return
            for k in range(min(n,mx),0,-1):
                for rest in gen(n-k,k):
                    yield (k,)+rest
        return list(gen(n,n))
    ok=True
    for n in range(0,17):
        for lam in parts(n):
            core,quot = core_quotient_4(lam)
            qsz = sum(sum(q) for q in quot)
            if sum(core) + 4*qsz != n:
                print("SIZE FAIL", lam, core, quot); ok=False
    print("core/quotient size self-check:", "OK" if ok else "FAIL")
selfcheck()

# ================= full residual table =================
from functools import lru_cache
@lru_cache(maxsize=None)
def Gtab(m):
    return psi_power_schur(m)

def G_of(part):
    """G_part(i)=(Re,Im) for part |- even. part may be ()."""
    n = sum(part)
    assert n % 2 == 0
    if n == 0: return (1,0)
    d = Gtab(n//2)
    return d.get(tuple(part), (0,0))

def run(MMAX=10):
    rows = []
    for m in range(1, MMAX+1):
        d = Gtab(m)
        for lam,(re,im) in d.items():
            vp = vpi(re,im)
            fl = hook_f(lam); v2f = v2(fl)
            core,quot = core_quotient_4(lam)
            gc = G_of(core); vpc = vpi(*gc)
            fc = hook_f(core); v2fc = v2(fc)
            rows.append(dict(m=m,lam=lam,re=re,im=im,vp=vp,v2f=v2f,
                             core=core,quot=quot,vpc=vpc,v2fc=v2fc))
    return rows

rows = run(10)
print("total shapes:", len(rows))

# Confirm unique full vanisher
van = [r['lam'] for r in rows if r['vp']==float('inf')]
print("full vanishers G_lam=0 (n<=20):", van)

# core=(2,2) family: is v_pi finite iff quotient nonempty?
print("\n-- core=(2,2) family --")
c22 = [r for r in rows if r['core']==(2,2)]
bad=[]
for r in c22:
    qsz = sum(sum(q) for q in r['quot'])
    fin = r['vp']!=float('inf')
    if (qsz>0) != fin: bad.append((r['lam'],r['quot'],r['vp']))
print("core=(2,2) count:", len(c22), " violations of (finite<=>nonempty quot):", bad)

# residual, EXCLUDING infinite cases (core=(2,2) gives vpc=inf)
def residual(r):
    if r['vpc']==float('inf') or r['vp']==float('inf'): return None
    return r['vp'] - r['v2f'] - r['vpc'] + r['v2fc']

# Group by ORDERED quotient and by UNORDERED multiset
from collections import defaultdict
def group_test(keyfn, label):
    groups = defaultdict(list)
    for r in rows:
        res = residual(r)
        if res is None: continue
        groups[keyfn(r)].append((r['lam'], res))
    nonconst=[]
    for k,v in groups.items():
        resset = set(x[1] for x in v)
        if len(resset)>1:
            nonconst.append((k, v[:6], sorted(resset)))
    print(f"\n[{label}] #groups={len(groups)}  #non-constant groups={len(nonconst)}")
    for k,sample,rs in nonconst[:12]:
        print("   quotient",k," residuals",rs," e.g.",sample[:4])
    return nonconst

nc_ord = group_test(lambda r: r['quot'], "ORDERED 4-quotient")
nc_un  = group_test(lambda r: tuple(sorted([q for q in r['quot']], key=lambda p:(sum(p),p))), "UNORDERED 4-quotient")
