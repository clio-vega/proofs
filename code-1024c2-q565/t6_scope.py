"""T6. SCOPE PREDICTED BY THE PROOF.

The proof of the core identity uses the sign eta = prod eps only through the fact that
eta is INVARIANT under the double flip (eps_a,eps_b) -> (-eps_a,-eps_b).  Any power of
eta therefore works.  Shifting theta by a constant c multiplies A(T) by prod_{x in T} x^c,
i.e. multiplies the term by eta^c.  So the theorem should hold for every

        theta = (c, c+1, ..., c+r-1),   c in Z

i.e. every theta whose rho-shifts m_j = theta_j+j-1 form an arithmetic progression of
common difference 2.  PRE-REGISTERED PREDICTION, tested here.

Negative control in the same run: theta NOT of this form should NOT vanish.
"""
import random
from fractions import Fraction as F
from eng import psi_pair, psi_head, admissible

random.seed(77)
def draw(lam, f, tries=500):
    for _ in range(tries):
        v = [F(random.randint(2,40), random.randint(1,9)) for _ in range(lam)]
        u = [F(random.choice([-1,1])*random.randint(2,40), random.randint(1,9)) for _ in range(f)]
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

print("A. shifted theta = (c,c+1,...,c+r-1)   [prediction: Psi = 0 whenever r >= f+2]")
ok = bad = 0
for lam in range(2,5):
    for f in range(0,3):
        for r in range(max(2,f+2), min(2*0+lam+f,6)+1):
            if r > lam+f: continue
            for c in [-3,-2,-1,0,1,2,3]:
                M = draw(lam,f)
                if M is None: continue
                theta = [c+j for j in range(r)]
                try: a = psi_pair(M, theta)
                except ZeroDivisionError: continue
                if a == 0: ok += 1
                else:
                    bad += 1
                    if bad <= 3: print("   FAIL", lam, f, r, c, a)
print(f"   zero in {ok}/{ok+bad} cells  (lam 2..4, f 0..2, r >= f+2, c in -3..3)")

print()
print("B. NEGATIVE CONTROL: theta NOT an AP of difference 2 in m  [must be NONZERO]")
ok2 = bad2 = 0
for lam in range(2,5):
    for f in range(0,3):
        for r in range(max(2,f+2), min(lam+f,5)+1):
            for _ in range(6):
                M = draw(lam,f)
                if M is None: continue
                while True:
                    theta = sorted(random.sample(range(-3,6), r))
                    m = [theta[j]+j for j in range(r)]
                    if len(set(m)) == r and any(m[j+1]-m[j] != 2 for j in range(r-1)): break
                try: a = psi_pair(M, theta)
                except ZeroDivisionError: continue
                if a != 0: ok2 += 1
                else:
                    bad2 += 1
                    if bad2 <= 4: print("   control fired: theta=", theta, " m=", m, " lam,f,r=", lam,f,r)
print(f"   nonzero in {ok2}/{ok2+bad2} draws  (control should stay green: nonzero everywhere)")
