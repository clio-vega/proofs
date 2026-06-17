"""
generalc_master.py -- prove session 2026-06-17

Verify the DERIVED general-c master valuation:

   R_i = G_{b+i}/G_0 ,   Delta(b+i) = (b+i) + 2 v2(R_i),
   v2(R_i) = v2(N_i) - c_i - v2(a-c+2) - [i>=2] v2(a-b+1) + v2(Pi_i) - E,

where (a=2P+b+c, m=P+b+c, w=b+c):
   c_i := v2(const_i) - v2(c!),   const_i from M_{b+i}=(b-c+1)*prod_{t=2}^i(b+t)*N_i/const_i (i>=2),
                                  M_{b+1}=(b-c+1)*(a-b+1)*N_1/const_1
   ell = ceil((w+3)/2),  E = (w-1)//2 if w odd else (w-2)//2,
   Pi_i = prod_{s=P+c-i+1}^{P+ell-1} s   (consecutive run, length L_i=ell-c+i-1).

Strategy: build N_i, const_i from symbolic fit; then check the master formula
against the EXACT R_i computed from the MN engine over a numeric range.
"""
import sympy as sp
from math import comb, factorial
from fractions import Fraction
from mn import Mj

a, b = sp.symbols('a b')

def collect(c, i, amax=64, bmax=28):
    pts = []
    for bb in range(c, bmax):
        for aa in range(bb, amax):
            if (aa + bb + c) % 2 != 0:
                continue
            m = (aa + bb + c) // 2
            j = bb + i
            if j > m:
                continue
            pts.append((aa, bb, Mj((aa, bb, c), j, m)))
    return pts

def fit_poly(pts, maxdeg):
    monos = [(p, q) for p in range(maxdeg + 1) for q in range(maxdeg + 1) if p + q <= maxdeg]
    if len(pts) < len(monos) + 3:
        return None
    A = sp.Matrix([[aa**p * bb**q for (p, q) in monos] for (aa, bb, _) in pts])
    y = sp.Matrix([M for (_, _, M) in pts])
    sol = A.solve_least_squares(y)
    if any(A * sol - y):
        return None
    expr = sum(sol[k] * a**p * b**q for k, (p, q) in enumerate(monos))
    return sp.factor(sp.expand(expr))

def extract(c, i):
    """return (N_i, const_i) with sign normalized so leading a-coeff>0; const_i rational."""
    M = fit_poly(collect(c, i), maxdeg=2 * c + 1)
    if M is None:
        return None, None, None
    bprod = (b - c + 1) * sp.prod([b + t for t in range(2, i + 1)])
    if i == 1:
        bprod = (b - c + 1) * (a - b + 1)
    # N_i*const_i_inv = M / bprod  ; we want M = bprod*N_i/const_i
    ratio = sp.cancel(M / bprod)   # = N_i/const_i
    # ratio is poly(a,b)/const. Extract integer content as const.
    num, den = sp.fraction(sp.together(ratio))
    num = sp.expand(num);
    poly = sp.Poly(num, a, b)
    coeffs = [sp.Rational(cc) for cc in poly.coeffs()]
    from sympy import gcd as sgcd
    g = 0
    for cc in coeffs:
        g = sp.igcd(g, sp.Integer(cc)) if cc==int(cc) else 1
    # const_i = den / g   (so N_i = num/g integer-primitive)
    const_i = sp.nsimplify(den / g)
    N_i = sp.expand(num / g)
    return N_i, const_i, M

def v2int(n):
    n = abs(int(n))
    if n == 0: return 999
    k = 0
    while n % 2 == 0: n//=2; k+=1
    return k

def v2run(lo, hi):
    """v2 of prod_{s=lo}^{hi} s ; empty if lo>hi -> 0."""
    tot = 0
    for s in range(lo, hi+1):
        tot += v2int(s)
    return tot

for c in [4, 5, 6]:
    print("="*70)
    print("c = %d" % c)
    print("="*70)
    info = {}
    for i in range(1, c+1):
        N_i, const_i, M = extract(c, i)
        if N_i is None:
            print("  i=%d: NO FIT"%i); continue
        c_i = v2int(int(factorial(c)/sp.Integer(const_i)*const_i)) # placeholder
        # c_i = v2(const_i) - v2(c!)
        ci = v2int(int(const_i)) - v2int(factorial(c)) if const_i==int(const_i) else None
        info[i] = (N_i, const_i, ci)
        print("  i=%d (k=%d): const_i=%s  c_i=v2(const)-v2(c!)=%s  N_i a-deg=%d"
              % (i, c-i, const_i, ci, sp.Poly(sp.expand(N_i),a).degree() if N_i.has(a) else 0))
    # --- verify master formula vs exact MN ---
    print("  --- verifying master v2(R_i) formula vs exact MN ---")
    w_parity_done=set()
    bad=0; checked=0
    for m in range(c+5, 40):
        for bb in range(2*c, m):
            aa = 2*m - c - bb
            if aa < bb: continue
            P = m - bb - c
            if P < 0: continue
            lam=(aa,bb,c)
            M0 = Mj(lam,0,m)
            if M0==0: continue
            w = bb+c
            ell = (w+3+1)//2 if (w+3)%2 else (w+3)//2   # ceil((w+3)/2)
            E = (w-1)//2 if w%2==1 else (w-2)//2
            for i in range(1,c+1):
                j=bb+i
                if j>m: continue
                try: Mbi = Mj(lam,j,m)
                except AssertionError: continue
                if Mbi==0: continue
                Ri = Fraction(comb(m,j))*Fraction(Mbi)/Fraction(M0)
                # exact v2 of Ri
                num=abs(Ri.numerator); den=Ri.denominator
                v2Ri = v2int(num)-v2int(den)
                # master prediction
                N_i, const_i, ci = info[i]
                Nval = int(N_i.subs({a:aa,b:bb}))
                vN = v2int(Nval)
                v_acm2 = v2int(2*P+bb+2)   # a-c+2
                v_abp1 = v2int(2*P+c+1)    # a-b+1
                # Pi_i = prod_{s=P+c-i+1}^{P+ell-1}
                vPi = v2run(P+c-i+1, P+ell-1)
                pred = vN - ci - v_acm2 - (v_abp1 if i>=2 else 0) + vPi - E
                checked+=1
                if pred != v2Ri:
                    bad+=1
                    if bad<=6:
                        print("     MISMATCH (a,b,m,i)=(%d,%d,%d,%d): pred=%d actual=%d"%(aa,bb,m,i,pred,v2Ri))
    print("  master formula: checked=%d  mismatches=%d"%(checked,bad))
