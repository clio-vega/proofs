"""
c5_certify.py -- prove session 2026-06-17

Close the c=5 boundary lemma concretely.  For each shape (a,b,5) in the box-interior
range, compute:
  (1) the TRUE Delta(b+i) via the MN engine,
  (2) a HAND LOWER BOUND that mirrors a clean proof:
        Delta(b+i) = (eta-k) - 2 v2(k!) + 2 v2(N_i) - 2 v2(a-c+2)
                     - 2 [k<=c-2] v2(a-b+1) + 2 v2(Pi_i)
      Even deficits (a-c+2 if b even; a-b+1 if i>=2) each cost the "+1" = -2.
      Their factor f (=(deficit)/2): if f in Pi_i range, absorb via Lemma P
      (remaining product v2 >= floor(L_i/2)-d_in); if out of range, compensate
      via v2(N_i) - v2(f) >= m_comp (a proven/grid surplus, >=0).
Goal: hand bound is (a) <= true Delta everywhere, (b) >= -theta+1 (strict).
theta = 0 (a even <-> b odd) ; theta = 3 (a odd <-> b even).
"""
from math import comb, factorial
import sympy as sp
from mn import Mj

a, b = sp.symbols('a b')
c = 5

def v2(n):
    n=abs(int(n))
    if n==0: return 999
    k=0
    while n%2==0: n//=2;k+=1
    return k

# --- explicit N_i^{(5)} (from fits, verified) ---
N = {
 5: sp.Integer(1),
 4: a*b+5*a-b**2-4*b-15,
 3: a**2*b**2+9*a**2*b+20*a**2-2*a*b**3-13*a*b**2-45*a*b-70*a+b**4+4*b**3+29*b**2+56*b+30,
 2: a**2*b**3+12*a**2*b**2+47*a**2*b+60*a**2-2*a*b**4-13*a*b**3-52*a*b**2-143*a*b-150*a+b**5+b**4+35*b**3+53*b**2-90,
 1: a**3*b**4+14*a**3*b**3+71*a**3*b**2+154*a**3*b+120*a**3-3*a**2*b**5-21*a**2*b**4-59*a**2*b**3-51*a**2*b**2+134*a**2*b+240*a**2+3*a*b**6+51*a*b**4+24*a*b**3-462*a*b**2-1296*a*b-1080*a-b**7+7*b**6-63*b**5+209*b**4-8*b**3-288*b**2-1296*b-2160,
}

def Pirange(P,bb,i):
    k=c-i; w=bb+c
    ell = -(-(w+3)//2)
    return (P+k+1, P+ell-1)   # lower, upper

def v2_run(lo,hi):
    return sum(v2(s) for s in range(lo,hi+1)) if hi>=lo else 0

def lemmaP_floor(L, d):
    # product of L consecutive ints has v2 >= floor(L/2); remove d factors -> >= floor(L/2)-d
    return L//2 - d

def v2_Ni(P, bb, i):
    """exact v2(N_i) from MN value of M_{b+i}:
       i>=2: M*c!*k! = (b-c+1) prod_{t=2}^i(b+t) N_i
       i=1 : M*c!*k! = (b-c+1)(a-b+1) N_1"""
    aa=2*P+bb+c; m=P+bb+c; k=c-i
    Mbi = Mj((aa,bb,c), bb+i, m)
    base = v2(Mbi) + v2(factorial(c)) + v2(factorial(k)) - v2(bb-c+1)
    if i==1:
        return base - v2(aa-bb+1)
    return base - sum(v2(bb+t) for t in range(2,i+1))

def hand_bound(P, bb, i):
    aa = 2*P+bb+c
    k = c-i; w=bb+c
    eta = 1 if w%2==1 else 2
    ell = -(-(w+3)//2)
    L = ell-k-1
    base = (eta-k) - 2*v2(factorial(k))
    lo,hi = Pirange(P,bb,i)
    # collect even deficits
    d_in = 0; out_facs = []
    # a-c+2 = 2P+b+2
    if bb%2==0:
        f1 = P + bb//2 + 1   # (a-c+2)/2
        base -= 2  # the +1
        if lo<=f1<=hi: d_in += 1
        else: out_facs.append(f1)
    # a-b+1 = 2P+6 (c=5 odd), present i>=2 (k<=3)
    if i>=2:
        f2 = P + 3
        base -= 2
        if lo<=f2<=hi: d_in += 1
        else: out_facs.append(f2)
    # Pi_i contribution after removing in-range deficit factors:
    piterm = 2*lemmaP_floor(L, d_in)
    # N_i compensation for out-of-range deficits + free surplus:
    comp = v2_Ni(P,bb,i) - sum(v2(f) for f in out_facs)   # must be >=0 (compensation)
    Nterm = 2*comp
    return base + piterm + Nterm, comp

def true_delta(P, bb, i):
    aa=2*P+bb+c; m=P+bb+c
    lam=(aa,bb,c)
    M0=Mj(lam,0,m);
    if M0==0: return None
    j=bb+i
    if j>m: return None
    try: Mj_=Mj(lam,j,m)
    except AssertionError: return None
    if Mj_==0: return None
    v0 = 0 + 2*v2(M0)              # val(0)=0+2 v2(C(m,0)M_0)=2 v2(M0)
    vj = j + 2*v2(comb(m,j)*Mj_)
    return vj - v0

# --- sweep ---
print("c=5 boundary certification (hand bound vs true Delta)")
print("b>=b0; b0 chosen so box interior fits. Using b>=12.")
results = {}  # (parity, i) -> [min_true_margin, min_hand_margin, min_comp]
bad_valid = 0; bad_suff = 0; checked=0
for P in range(0, 60):
    for bb in range(12, 90):
        aa=2*P+bb+c
        if aa<bb: continue
        m=P+bb+c
        if m>140: continue
        par = 'aEVEN' if aa%2==0 else 'aODD'
        theta = 0 if aa%2==0 else 3
        for i in range(1,c+1):
            td = true_delta(P,bb,i)
            if td is None: continue
            hb, comp = hand_bound(P,bb,i)
            checked+=1
            key=(par,i)
            cur = results.get(key, [999,999,999])
            cur[0]=min(cur[0], td+theta)   # true margin above floor
            cur[1]=min(cur[1], hb+theta)
            cur[2]=min(cur[2], comp)
            results[key]=cur
            if hb > td:        # bound must be VALID (<= true)
                bad_valid+=1
                if bad_valid<=5: print("  !! bound>true (P=%d,b=%d,i=%d): hb=%d td=%d"%(P,bb,i,hb,td))
            if not (hb > -theta):   # bound must be SUFFICIENT (strict)
                bad_suff+=1
                if bad_suff<=5: print("  !! bound not strict (P=%d,b=%d,i=%d): hb=%d theta=%d"%(P,bb,i,hb,theta))

print("checked=%d  bound-invalid(hb>true)=%d  bound-insufficient(hb<=-theta)=%d"%(checked,bad_valid,bad_suff))
print("\nPer (parity,index): min(true+theta), min(hand+theta), min(comp=v2N-out-deficits)")
for par in ['aEVEN','aODD']:
    for i in range(c,0,-1):
        if (par,i) in results:
            r=results[(par,i)]
            print("  %s i=%d (k=%d): true margin>=%d, hand margin>=%d, comp>=%d"%(par,i,c-i,r[0],r[1],r[2]))
