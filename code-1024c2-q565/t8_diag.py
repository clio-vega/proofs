"""T8. Diagnose the 21 firings of the mixed-parity control.

Two mechanisms are already known to be theta-INDEPENDENT:
  (G) global negation: if f=0 then M=-M, and B(-T)=B(T) while A(-T)=(-1)^{C(r,2)}A(T)
      for EVERY theta, so Psi = (-1)^{C(r,2)} Psi, hence Psi=0 when r = 2,3 mod 4.
Ask: are the firings all of type (G)?  And are the non-(G) ones IDENTICAL on the
stratum (several independent draws) or accidental (one draw)?
"""
import random
from fractions import Fraction as F
from eng import psi_pair, admissible

random.seed(1234)
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

def theta_from_m(m): return [m[j]-j for j in range(len(m))]

# (G) check: f=0 and C(r,2) odd, ARBITRARY theta
print("(G)  f=0 and r(r-1)/2 odd, arbitrary mixed-parity m  [predict Psi=0 for every theta]")
ok=bad=0
for lam in range(2,5):
    for r in [2,3,6]:
        if r > lam or (r*(r-1)//2) % 2 == 0: continue
        for _ in range(12):
            M = draw(lam,0)
            m = sorted(random.sample(range(-5,9), r))
            try: a = psi_pair(M, theta_from_m(m))
            except ZeroDivisionError: continue
            ok += (a==0); bad += (a!=0)
print(f"     zero in {ok}/{ok+bad}")

print()
print("(G') f=0 and r(r-1)/2 EVEN, mixed-parity m  [predict NONZERO: no mechanism]")
ok=bad=0
for lam in range(4,6):
    for r in [4,5]:
        if r > lam: continue
        for _ in range(12):
            M = draw(lam,0)
            if M is None: continue
            while True:
                m = sorted(random.sample(range(-5,9), r))
                if len(set(x%2 for x in m))==2: break
            try: a = psi_pair(M, theta_from_m(m))
            except ZeroDivisionError: continue
            ok += (a!=0); bad += (a==0)
print(f"     nonzero in {ok}/{ok+bad}")

print()
print("non-(G) firing: lam=3, f=1, r=4, m=[-1,1,4,7]  -- identical or accidental?")
m = [-1,1,4,7]
vals = []
for _ in range(10):
    M = draw(3,1)
    if M is None: continue
    try: a = psi_pair(M, theta_from_m(m))
    except ZeroDivisionError: continue
    vals.append(a)
print(f"     {sum(1 for x in vals if x==0)}/{len(vals)} draws give zero")
print(f"     sample values: {[str(x)[:28] for x in vals[:4]]}")
