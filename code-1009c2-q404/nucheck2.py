import sympy as sp, mpmath as mp
t=sp.Symbol('t'); mp.mp.dps=40
def nu_and_M(p):
    cs=[mp.mpf(int(c)) for c in sp.Poly(p,t).all_coeffs()]
    rts=mp.polyroots(cs,maxsteps=600,extraprec=800) if len(cs)>1 else []
    nz=[r for r in rts if abs(r)>mp.mpf('1e-30')]
    nv=mp.fprod([max(abs(r),1/abs(r)) for r in nz]) if nz else mp.mpf(1)
    M =abs(cs[0])*mp.fprod([max(mp.mpf(1),abs(r)) for r in rts]) if rts else abs(cs[0])
    return nv,M
bad=0; mins=[]
for b in range(2,51):
    nv,M=nu_and_M(t**b-t**(b-1)+2)
    if abs(M-2)>mp.mpf('1e-25') or abs(nv-2)>mp.mpf('1e-25'):
        bad+=1; print('  FAIL E_%d M=%s nu=%s'%(b,mp.nstr(M,12),mp.nstr(nv,12)))
    cs=[mp.mpf(int(c)) for c in sp.Poly(t**b-t**(b-1)+2,t).all_coeffs()]
    mins.append(min(abs(r) for r in mp.polyroots(cs,maxsteps=600,extraprec=800)))
print('M(E_b) = nu(E_b) = 2 exactly, 2<=b<=50          : %d failures'%bad)
print('min root modulus of E_b >= 1 (no root in the open disc), 2<=b<=50 : %d failures'
      %sum(1 for m in mins if m < 1-mp.mpf('1e-25')))
print('  smallest min|root| observed =', mp.nstr(min(mins),12))
# reciprocality, exact
bad=0; rows=[]
for b in range(1,61):
    Db=t**b-t**(b-1)+1; Dst=sp.expand(t**b*Db.subs(t,1/t))
    d1=sp.expand(Dst-Db); d2=sp.expand(Dst+Db)
    if b>=3 and (d1==0 or d2==0): bad+=1
    if b==2 and d1!=0: bad+=1
    if b<=4: rows.append('  b=%d: D_b^* = %s ;  D_b^* - D_b = %s'%(b,Dst,d1))
print()
print('D_b^* != +-D_b for 3<=b<=60, and D_2^* = D_2   : %d failures'%bad)
for r in rows: print(r)
print('  closed form claimed in the proof: D_b^* - D_b = t^{b-1} - t ;  check b=1..60:',
      sum(1 for b in range(1,61)
          if sp.expand(sp.expand(t**b*(t**b-t**(b-1)+1).subs(t,1/t))-(t**b-t**(b-1)+1))
             != sp.expand(t**(b-1)-t)), 'failures')
# E_b reciprocality for completeness
print('  E_b^* - E_b for b=2,3:', [sp.expand(sp.expand(t**b*(t**b-t**(b-1)+2).subs(t,1/t))-(t**b-t**(b-1)+2)) for b in (2,3)])
