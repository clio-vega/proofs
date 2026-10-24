"""T11.  The four "disagreements" of T10 diagnosed, and the theorem sharpened.

All four are L = s+1 with nu odd.  Reason: substituting x_a -> -x_a and reindexing
eps_a -> -eps_a gives Z(-x_a) = (-1)^nu Z.  So for nu ODD, Z is odd in each x_a,
hence x_a | P for every a, hence  Delta * prod_a x_a  divides P.  Its degree in x_1 is
2(L-1)+1, and the Step-1 bound is still N-1 = L+s-1, so

    2(L-1)+1 <= L+s-1   =>   L <= s,      i.e.   L >= s+1  forces  P = 0.

Since nu = N-1+delta and L=s+1 gives N=2s+1, nu = 2s+delta = delta (mod 2).  So:

    SHARPENED THEOREM:  all m_j of parity delta  =>  Psi = 0 whenever r >= f + 2 - delta.

i.e. r >= f+2 for all-EVEN m, and r >= f+1 for all-ODD m.  TESTED HERE ON Psi ITSELF.
"""
import random
from fractions import Fraction as F
from eng import psi_pair, psi_head, admissible

random.seed(5150)
def draw(lam, f, tries=900):
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
def odd_m(r):
    while True:
        k = sorted(random.sample(range(-3,7), r))
        m = [2*x+1 for x in k]
        if len(set(m)) == r: return m
def even_m(r):
    while True:
        k = sorted(random.sample(range(-3,7), r))
        m = [2*x for x in k]
        if len(set(m)) == r: return m

print("A.  all m_j ODD at the NEW cell r = f+1  [prediction: Psi = 0]")
ok=bad=0; fails=[]
for lam in range(1,6):
    for f in range(0,4):
        r = f+1
        if r > lam+f or r < 1: continue
        for _ in range(8):
            M = draw(lam,f)
            if M is None: continue
            m = odd_m(r)
            try:
                a = psi_pair(M, th(m)); b = psi_head(M, th(m))
            except ZeroDivisionError: continue
            assert a == b, "engines disagree"
            if a == 0: ok += 1
            else:
                bad += 1
                if len(fails)<3: fails.append((lam,f,r,m,a))
print(f"    zero in {ok}/{ok+bad}   (lam<=5, f<=3, r=f+1, both engines agree on every draw)")
for x in fails: print("    FAIL", x)

print()
print("B.  CONTROL: all m_j EVEN at r = f+1  [prediction: NONZERO -- delta=0 needs r>=f+2]")
ok=bad=0; fired=[]
for lam in range(1,6):
    for f in range(0,4):
        r = f+1
        if r > lam+f or r < 1: continue
        for _ in range(8):
            M = draw(lam,f)
            if M is None: continue
            m = even_m(r)
            try: a = psi_pair(M, th(m))
            except ZeroDivisionError: continue
            if a != 0: ok += 1
            else:
                bad += 1
                if len(fired)<4: fired.append((lam,f,r,m))
print(f"    nonzero in {ok}/{ok+bad}")
for x in fired: print("    control fired:", x)

print()
print("C.  CONTROL: all m_j ODD at r = f  (one BELOW the new threshold) [prediction: NONZERO]")
ok=bad=0; fired=[]
for lam in range(1,6):
    for f in range(1,4):
        r = f
        if r > lam+f: continue
        for _ in range(8):
            M = draw(lam,f)
            if M is None: continue
            m = odd_m(r)
            try: a = psi_pair(M, th(m))
            except ZeroDivisionError: continue
            if a != 0: ok += 1
            else:
                bad += 1
                if len(fired)<4: fired.append((lam,f,r,m))
print(f"    nonzero in {ok}/{ok+bad}")
for x in fired: print("    control fired:", x)

print()
print("D.  the full sharpened criterion r >= f+2-delta, swept")
ok=bad=0
for lam in range(1,6):
    for f in range(0,4):
        for r in range(1, min(lam+f,6)+1):
            for delta in (0,1):
                if r < f+2-delta: continue
                for _ in range(4):
                    M = draw(lam,f)
                    if M is None: continue
                    m = odd_m(r) if delta else even_m(r)
                    try: a = psi_pair(M, th(m))
                    except ZeroDivisionError: continue
                    ok += (a==0); bad += (a!=0)
print(f"    zero in {ok}/{ok+bad}   (lam<=5, f<=3, r<=6, both delta)")
