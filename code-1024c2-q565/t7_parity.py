"""T7. The control firing of T6 diagnosed.

If all m_j = theta_j+j-1 are congruent mod 2, say m_j = 2k_j+delta, then
   det(x_i^{m_j}) = (prod x_i)^delta * det((x_i^2)^{k_j}),
so    A(T) = (prod_{x in T} x)^delta * s_mu(x_1^2,...,x_r^2) * prod_{i<k}(x_i+x_k).
Within a class the middle factor is eps-INVARIANT (a_i^2 = v_i^2) and the first is
eta^delta.  So the proof applies verbatim with eta^{N-1} -> eta^{N-1+delta}.

SHARPENED SCOPE:  m_j pairwise distinct and ALL OF THE SAME PARITY.
Sharpened control: m of MIXED parity must be NONZERO.
"""
import random, itertools
from fractions import Fraction as F
from eng import psi_pair, admissible, alternant

random.seed(99)
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

def theta_from_m(m):   # m_j = theta_j + j - 1
    return [m[j]-j for j in range(len(m))]

print("A. A(T) = (prod x)^delta * s_mu(x^2) * prod(x_i+x_k)  when all m_j = delta mod 2")
bad = ok = 0
for r in range(1,5):
    for _ in range(60):
        delta = random.choice([0,1])
        k = sorted(random.sample(range(0,7), r))
        m = [2*kk+delta for kk in k]
        x = [F(random.randint(-20,20), random.randint(1,5)) for _ in range(r)]
        if len(set(x)) != r or any(xx==0 for xx in x) or len(set(abs(xx) for xx in x)) != r: continue
        A = alternant(x, m)
        pre = F(1)
        for xx in x: pre *= xx**delta
        sq = alternant([xx**2 for xx in x], k)     # s_mu in the squares
        prod = F(1)
        for i in range(r):
            for j in range(i+1,r): prod *= (x[i]+x[j])
        if A == pre*sq*prod: ok += 1
        else:
            bad += 1
            if bad<=2: print("   FAIL", m, x, A, pre*sq*prod)
print(f"   {ok}/{ok+bad} exact  (r<=4, random x with distinct |x|)")

print()
print("B. SAME PARITY m, arbitrary spacing  [prediction: Psi=0 when r >= f+2]")
ok = bad = 0
for lam in range(2,5):
    for f in range(0,3):
        for r in range(max(2,f+2), min(lam+f,5)+1):
            for _ in range(10):
                M = draw(lam,f)
                if M is None: continue
                delta = random.choice([0,1])
                m = [2*kk+delta for kk in sorted(random.sample(range(-3,7), r))]
                try: a = psi_pair(M, theta_from_m(m))
                except ZeroDivisionError: continue
                if a == 0: ok += 1
                else:
                    bad += 1
                    if bad<=3: print("   FAIL  m=",m," lam,f,r=",lam,f,r, " Psi=",a)
print(f"   zero in {ok}/{ok+bad} draws")

print()
print("C. SHARPENED CONTROL: m of MIXED parity  [must be NONZERO]")
ok = bad = 0; firings = []
for lam in range(2,5):
    for f in range(0,3):
        for r in range(max(2,f+2), min(lam+f,5)+1):
            for _ in range(10):
                M = draw(lam,f)
                if M is None: continue
                while True:
                    m = sorted(random.sample(range(-4,8), r))
                    if len(set(mm % 2 for mm in m)) == 2: break
                try: a = psi_pair(M, theta_from_m(m))
                except ZeroDivisionError: continue
                if a != 0: ok += 1
                else:
                    bad += 1
                    firings.append((m, lam, f, r))
print(f"   nonzero in {ok}/{ok+bad} draws")
for fr in firings[:6]: print("   control fired: m=", fr[0], " lam,f,r=", fr[1:])
