"""E_b (the a=b diagonal) and the OFF-diagonal cells of Theorem C."""
import sympy as sp, mpmath as mp, sys
t=sp.Symbol('t'); mp.mp.dps=40

def E(b): return sp.Poly(t**b - t**(b-1) + 2, t)

print("--- E_b = t^b - t^{b-1} + 2 ---")
print(" b | E_b(0) | #roots on |z|=1 | those roots | min|a| | max|a| | nu(E_b) | M(E_b)")
badconst=0; badcirc=0; badnu=0; BM=30
for b in range(1,BM+1):
    p=E(b); cs=[mp.mpf(int(c)) for c in p.all_coeffs()]
    rts=mp.polyroots(cs,maxsteps=400,extraprec=600)
    onc=[r for r in rts if abs(abs(r)-1)<mp.mpf('1e-25')]
    nv=mp.fprod([max(abs(r),1/abs(r)) for r in rts])
    M=mp.fprod([max(mp.mpf(1),abs(r)) for r in rts])
    c0=p.as_expr().subs(t,0)
    # CLAIM: E_b(0) = 2 for b>=2 and = 1 for b=1
    if c0 != (2 if b>=2 else 1): badconst+=1
    # CLAIM: roots on |z|=1 are exactly {-1}, occurring iff b is odd
    want = 1 if b%2==1 else 0
    if len(onc)!=want: badcirc+=1
    if onc and any(abs(r+1)>mp.mpf('1e-20') for r in onc): badcirc+=1
    # CLAIM: nu(E_b) >= 2 for b>=2 ; nu(E_1)=1
    if b>=2 and not nv>=2-mp.mpf('1e-20'): badnu+=1
    if b==1 and abs(nv-1)>mp.mpf('1e-20'): badnu+=1
    if b<=10 or b==BM:
        print("%2d | %s | %d | %s | %s | %s | %s | %s"%(b,c0,len(onc),
            ",".join(mp.nstr(r,5) for r in onc) if onc else "-",
            mp.nstr(min(abs(r) for r in rts),6), mp.nstr(max(abs(r) for r in rts),6),
            mp.nstr(nv,10), mp.nstr(M,10)))
print()
print("E_b, b=1..%d:"%BM)
print("  E_b(0) = 2 for b>=2, = 1 for b=1 (exact)                 : %d failures"%badconst)
print("  unit-circle roots are exactly {-1}, and only for b odd   : %d failures"%badcirc)
print("  nu(E_b) >= 2 for b>=2 ; nu(E_1) = 1                      : %d failures"%badnu)
print("  E_1 =", sp.factor(E(1).as_expr()), " E_2 =", sp.factor(E(2).as_expr()), " E_3 =", sp.factor(E(3).as_expr()))

print()
print("--- off-diagonal cells of Theorem C: are they product forms? ---")
# middle case 1<=m<=b-1 :  (t-1)t^{b-1} + (t-1)t^{b-m-1}  ?=  -t^{b-m-1}(1-t)(1-t^{2m})/(1-t^m)
bad=0; tot=0
for b in range(1,16):
    for m in range(1,b):
        lhs = sp.expand((t-1)*t**(b-1) + (t-1)*t**(b-m-1))
        rhs = sp.expand(-t**(b-m-1)*(1-t)*sp.cancel((1-t**(2*m))/(1-t**m)))
        tot+=1
        if sp.simplify(lhs-rhs)!=0:
            bad+=1
            if bad<4: print("  FAIL b=%d m=%d  %s  vs  %s"%(b,m,lhs,rhs))
print("  middle case = -t^{b-m-1}(1-t)(1-t^{2m})/(1-t^m) : %d (b,m) pairs, %d failures"%(tot,bad))
# and the m>b case
bad2=0;tot2=0
for b in range(1,16):
    lhs=sp.expand((t-1)*t**(b-1)); rhs=sp.expand(-t**(b-1)*(1-t)); tot2+=1
    if sp.simplify(lhs-rhs)!=0: bad2+=1
print("  m>b case   = -t^{b-1}(1-t)                      : %d values, %d failures"%(tot2,bad2))
# all roots of the off-diagonal cells on {0} u S^1 ?
badu=0; totu=0; examples=[]
for b in range(1,16):
    for m in list(range(1,b))+[None]:
        f = (t-1)*t**(b-1) + ((t-1)*t**(b-m-1) if m is not None else 0)
        p=sp.Poly(sp.expand(f),t); totu+=1
        cs=[mp.mpf(int(c)) for c in p.all_coeffs()]
        rts=mp.polyroots(cs,maxsteps=400,extraprec=600) if p.degree()>0 else []
        off=[r for r in rts if abs(r)>mp.mpf('1e-25') and abs(abs(r)-1)>mp.mpf('1e-20')]
        if off:
            badu+=1; examples.append((b,m,[mp.nstr(r,6) for r in off]))
print("  every root lies in {0} u S^1 (unitary)          : %d cells, %d failures"%(totu,badu))
for e in examples[:3]: print("   ",e)
# POSITIVE CONTROL for the unitarity scan: feed it the diagonal cells, which must FAIL
badc=0; totc=0
for b in range(3,16):
    for f in [t**b-t**(b-1)+1, t**b-t**(b-1)+2]:
        p=sp.Poly(f,t); totc+=1
        cs=[mp.mpf(int(c)) for c in p.all_coeffs()]
        rts=mp.polyroots(cs,maxsteps=400,extraprec=600)
        off=[r for r in rts if abs(r)>mp.mpf('1e-25') and abs(abs(r)-1)>mp.mpf('1e-20')]
        if off: badc+=1
print("  CONTROL: same scan on the DIAGONAL cells D_b,E_b (b=3..15): %d cells, %d non-unitary  [must equal %d]"%(totc,badc,totc))
