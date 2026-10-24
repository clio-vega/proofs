"""T9. (G) corrected, and a scan of parity patterns.

A(-T) = det((-x_i)^{m_j}) / prod((-x_k)+x_i) = (-1)^{sum m_j} / (-1)^{C(r,2)} * A(T),
so for f=0 (M=-M), since B(-T)=B(T):   Psi = (-1)^{sum_j m_j + C(r,2)} Psi.
=> (G) Psi=0 whenever f=0 and sum_j m_j + r(r-1)/2 is ODD, for EVERY theta.
(I first wrote this with the sum_j m_j factor missing.  The control caught it.)
"""
import random, itertools
from fractions import Fraction as F
from eng import psi_pair, admissible

random.seed(31337)
def draw(lam, f, tries=800):
    for _ in range(tries):
        v = [F(random.randint(2,60), random.randint(1,11)) for _ in range(lam)]
        u = [F(random.choice([-1,1])*random.randint(2,60), random.randint(1,11)) for _ in range(f)]
        M = []
        for x in v: M += [x,-x]
        M += u
        if not admissible(M, require_recip_free=True): continue
        e = F(1)
        for x in M: e *= x
        if e in (1,-1): continue
        if sum(1 for a in range(len(M)) for b in range(a+1,len(M)) if M[a]==-M[b]) != lam: continue
        return M
    return None
def th(m): return [m[j]-j for j in range(len(m))]

print("(G) corrected: f=0, sum(m)+C(r,2) ODD, arbitrary m  [predict Psi=0]")
ok=bad=0
for lam in range(2,6):
    for r in range(2,min(lam,5)+1):
        for _ in range(14):
            M = draw(lam,0)
            if M is None: continue
            for _try in range(50):
                m = sorted(random.sample(range(-5,9), r))
                if (sum(m)+r*(r-1)//2) % 2 == 1: break
            else: continue
            try: a = psi_pair(M, th(m))
            except ZeroDivisionError: continue
            ok += (a==0); bad += (a!=0)
print(f"    zero in {ok}/{ok+bad}")

print()
print("(G) control: f=0, sum(m)+C(r,2) EVEN, m MIXED parity  [predict NONZERO]")
ok=bad=0
for lam in range(3,6):
    for r in range(2,min(lam,5)+1):
        for _ in range(14):
            M = draw(lam,0)
            if M is None: continue
            for _try in range(80):
                m = sorted(random.sample(range(-5,9), r))
                if (sum(m)+r*(r-1)//2) % 2 == 0 and len(set(x%2 for x in m))==2: break
            else: continue
            try: a = psi_pair(M, th(m))
            except ZeroDivisionError: continue
            ok += (a!=0); bad += (a==0)
print(f"    nonzero in {ok}/{ok+bad}")

print()
print("PARITY-PATTERN SCAN:  lam=3, f=1, r=4  (so r>=f+2).  For each mod-2 pattern of m,")
print("is Psi identically 0 on the stratum (all 6 independent draws)?")
res = {}
for pat in itertools.product([0,1], repeat=4):
    zs = 0; nd = 0
    for _ in range(6):
        M = draw(3,1)
        if M is None: continue
        for _try in range(200):
            m = sorted(random.sample(range(-6,10), 4))
            if tuple(x % 2 for x in m) == pat: break
        else: continue
        try: a = psi_pair(M, th(m))
        except ZeroDivisionError: continue
        nd += 1; zs += (a==0)
    if nd: res[pat] = (zs, nd)
for pat,(zs,nd) in sorted(res.items()):
    same = len(set(pat))==1
    print(f"    m mod 2 = {pat}  (uniform={same})   zero in {zs}/{nd}")
