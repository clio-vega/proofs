import sympy as sp, mpmath as mp
t=sp.Symbol('t'); mp.mp.dps=40
def nu_and_M(p):
    cs=[mp.mpf(int(c)) for c in sp.Poly(p,t).all_coeffs()]
    rts=mp.polyroots(cs,maxsteps=600,extraprec=800) if len(cs)>1 else []
    nz=[r for r in rts if abs(r)>mp.mpf('1e-30')]          # SKIP the zero roots
    nv=mp.fprod([max(abs(r),1/abs(r)) for r in nz]) if nz else mp.mpf(1)
    M =abs(cs[0])*mp.fprod([max(mp.mpf(1),abs(r)) for r in rts]) if rts else abs(cs[0])
    return nv,M,len(rts)-len(nz)
bad=0; Ms=[]
for b in range(1,61):
    nv,M,nzero=nu_and_M(t**b-t**(b-1)+1)
    # nu = M^2 holds because D_b is monic with |D_b(0)|=1 (for b>=2); D_1=t has a zero root
    if b>=2 and abs(nv-M**2)>mp.mpf('1e-25'): bad+=1; print('  FAIL nu=M^2 at b=%d'%b)
    Ms.append((b,M))
print('nu(D_b) = M(D_b)^2 for 2<=b<=60 : %d failures'%bad)
th0=mp.findroot(lambda x: x**3-x-1, 1.3)
print('theta_0 (plastic number, root of x^3-x-1) =', mp.nstr(th0,12))
m3=[M for b,M in Ms if b==3][0]
print('M(D_3) =', mp.nstr(m3,12), '   M(D_3) - theta_0 =', mp.nstr(m3-th0,4))
sub=[(b,M) for b,M in Ms if b>=3]
print('min over 3<=b<=60 : b=%d  M=%s'%min(sub,key=lambda z:z[1])[::-1][::-1][0:2] if False else '')
mn=min(sub,key=lambda z:z[1]); mx=max(sub,key=lambda z:z[1])
print('min M(D_b), 3<=b<=60 : b=%d, M=%s'%(mn[0],mp.nstr(mn[1],12)))
print('max M(D_b), 3<=b<=60 : b=%d, M=%s'%(mx[0],mp.nstr(mx[1],12)))
print('all M(D_b) >= theta_0 for b>=3 :', all(M>=th0-mp.mpf('1e-25') for b,M in sub))
print('M(D_b) for b=50..60:', [ (b,mp.nstr(M,10)) for b,M in Ms if b>=50 ])
# the brief's table column, recomputed: max single root modulus vs M
print()
print(' b | max single |root| | #roots off S^1 | M(D_b)')
for b in [3,5,9,15,25,40,60]:
    cs=[mp.mpf(int(c)) for c in sp.Poly(t**b-t**(b-1)+1,t).all_coeffs()]
    rts=mp.polyroots(cs,maxsteps=600,extraprec=800)
    off=[r for r in rts if abs(abs(r)-1)>mp.mpf('1e-20')]
    nv,M,_=nu_and_M(t**b-t**(b-1)+1)
    print('%2d | %s | %d | %s'%(b,mp.nstr(max(abs(r) for r in rts),8),len(off),mp.nstr(M,10)))
# E_b: M = nu = 2 exactly?
bad=0
for b in range(2,61):
    nv,M,_=nu_and_M(t**b-t**(b-1)+2)
    if abs(M-2)>mp.mpf('1e-25') or abs(nv-2)>mp.mpf('1e-25'): bad+=1; print('  FAIL E_%d M=%s nu=%s'%(b,mp.nstr(M,12),mp.nstr(nv,12)))
print()
print('M(E_b) = nu(E_b) = 2 exactly, 2<=b<=60 : %d failures'%bad)
# reciprocality: D_b^* - D_b and D_b^* + D_b
bad=0
for b in range(1,61):
    Db=t**b-t**(b-1)+1; Dst=sp.expand(t**b*Db.subs(t,1/t))
    d1=sp.expand(Dst-Db); d2=sp.expand(Dst+Db)
    if b>=3 and (d1==0 or d2==0): bad+=1
    if b==2 and d1!=0: bad+=1
    if b<=4: print('  b=%d: D_b^* = %s ;  D_b^*-D_b = %s'%(b,Dst,d1))
print('D_b^* != +-D_b for 3<=b<=60, and D_2^* = D_2 : %d failures'%bad)
