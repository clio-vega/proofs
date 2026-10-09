"""Two mechanisms for M(D_b), separated.  The Jensen contour integral has a log
singularity exactly when Phi6 | D_b (b = 2 mod 6), so deflate Phi6 EXACTLY first;
M(Phi6) = 1, so this changes no value.  Tolerance is justified by comparing two
quadrature settings against each other, not by loosening until green."""
import sympy as sp, mpmath as mp, sys
t = sp.Symbol('t'); mp.mp.dps = 40
Phi6 = sp.Poly(t**2 - t + 1, t)
def D(b): return sp.Poly(t**b - t**(b-1) + 1, t)

def M_roots(poly):
    cs=[mp.mpf(int(c)) for c in poly.all_coeffs()]
    if len(cs)==1: return abs(cs[0])
    rts=mp.polyroots(cs,maxsteps=400,extraprec=600)
    return abs(cs[0])*mp.fprod([max(mp.mpf(1),abs(r)) for r in rts])

def M_jensen(poly, maxdegree):
    cs=[int(c) for c in poly.all_coeffs()]
    def f(th):
        z=mp.e**(2j*mp.pi*th); acc=mp.mpc(0)
        for c in cs: acc=acc*z+c
        return mp.log(abs(acc))
    return mp.e**mp.quad(f,[0,1],maxdegree=maxdegree)

BMAX=int(sys.argv[1]) if len(sys.argv)>1 else 40
bad_val=0; bad_conv=0; bad_nuM=0; worst_val=mp.mpf(0); worst_conv=mp.mpf(0)
for b in range(1,BMAX+1):
    p=D(b)
    # exact deflation of the only possible cyclotomic factor
    if b%6==2:
        q,r=sp.div(p,Phi6,t); assert sp.expand(r.as_expr())==0, b
        pd=sp.Poly(q,t)
    else:
        pd=p
    Mr=M_roots(pd)
    J1=M_jensen(pd,8); J2=M_jensen(pd,12)
    conv=abs(J1-J2)                  # the quadrature's own convergence estimate
    err=abs(Mr-J2)
    worst_conv=max(worst_conv,conv); worst_val=max(worst_val,err)
    if conv>mp.mpf('1e-15'): bad_conv+=1
    if err>mp.mpf('1e-15'): bad_val+=1
    # nu = M^2 for monic f with |f(0)|=1 (D_b and its Phi6-deflation both qualify)
    if abs(M_roots(p)**2 - (lambda q: q)(None) if False else 0)>0: pass
print("M(D_b) two ways, b=1..%d, Phi6 deflated exactly where b=2 mod 6:"%BMAX)
print("  quadrature self-convergence  max|J(deg8)-J(deg12)| = %s   (%d of %d above 1e-15)"%(mp.nstr(worst_conv,4),bad_conv,BMAX))
print("  mechanism disagreement       max|M_roots - M_Jensen| = %s   (%d of %d above 1e-15)"%(mp.nstr(worst_val,4),bad_val,BMAX))
# POSITIVE CONTROL: feed the Jensen arm a polynomial with a DIFFERENT Mahler measure
ctrl=sp.Poly(t**3-t-1,t)   # the plastic number, M = 1.3247...
print("  control: M_roots(t^3-t-1) = %s   M_Jensen(t^3-t-1) = %s"%(mp.nstr(M_roots(ctrl),12),mp.nstr(M_jensen(ctrl,12),12)))
print("  control (planted wrong target): |M_roots(t^3-t-1) - M_Jensen(D_5)| = %s  <-- must be LARGE"%(
      mp.nstr(abs(M_roots(ctrl)-M_jensen(D(5),12)),6)))
# nu = M^2 check, stated separately
bad=0
for b in range(1,BMAX+1):
    p=D(b); cs=[mp.mpf(int(c)) for c in p.all_coeffs()]
    rts=mp.polyroots(cs,maxsteps=400,extraprec=600)
    nv=mp.fprod([max(abs(r),1/abs(r)) for r in rts])
    if abs(nv-M_roots(p)**2)>mp.mpf('1e-30'): bad+=1
print("  nu(D_b) = M(D_b)^2  (valid because D_b is monic with |D_b(0)|=1): %d of %d failures"%(bad,BMAX))
