"""EXACT test of 'monomial times a product of cyclotomics', by rational factorisation.
   Independent of all root-modulus numerics."""
import sympy as sp, sys
from sympy.polys.specialpolys import cyclotomic_poly
t=sp.Symbol('t')

_CYC={}
def is_cyclotomic(p):
    """p a monic irreducible in Z[t]: is p = Phi_e for some e?"""
    d=p.degree()
    if d not in _CYC:
        # all e with totient(e) = d; totient(e) >= sqrt(e/2) so e <= 2 d^2 + 2 is a safe cap
        _CYC[d]=[sp.Poly(cyclotomic_poly(e,t),t) for e in range(1,2*d*d+3) if sp.totient(e)==d]
    return any(p==c for c in _CYC[d])

def unitary_exact(expr):
    """Return (is_monomial_times_cyclotomics, [non-cyclotomic irreducible factors])."""
    p=sp.Poly(sp.expand(expr),t)
    unit,facs=sp.factor_list(p.as_expr(), t)
    bad=[]
    for f,mult in facs:
        pf=sp.Poly(f,t)
        if pf.as_expr()==t: continue
        lc=pf.all_coeffs()[0]
        if lc!=1: pf=sp.Poly(sp.expand(pf.as_expr()/lc),t)
        if not is_cyclotomic(pf): bad.append((sp.factor(pf.as_expr()),mult))
    return (len(bad)==0), bad

if __name__=='__main__':
    # self-test of the instrument, both arms
    print("instrument self-test:")
    for e,lab in [(t,'t'),(t**2-t+1,'Phi_6'),(t**6-1,'t^6-1'),(1-t**12,'1-t^12'),
                  ((t-1)*t**4*(t**3+1),'(t-1)t^4(t^3+1)'),
                  (t**3-t-1,'t^3-t-1 (plastic, NOT cyclotomic)'),
                  (t**2-t+2,'t^2-t+2 (= E_2, NOT cyclotomic)'),
                  (2*t,'2t (monomial with a rational scalar)')]:
        ok,bad=unitary_exact(e)
        print("   %-34s -> %s   %s"%(lab,"CYCLOTOMIC x MONOMIAL" if ok else "NOT",bad if bad else ""))
    print()
    BM=int(sys.argv[1]) if len(sys.argv)>1 else 40
    # D_b
    nd=0; cycpart=[]
    for b in range(1,BM+1):
        ok,bad=unitary_exact(t**b-t**(b-1)+1)
        expect = (b<=2)
        if ok!=expect: nd+=1; print("  MISMATCH D_%d ok=%s expect=%s"%(b,ok,expect))
        unit,facs=sp.factor_list(sp.expand(t**b-t**(b-1)+1),t)
        cyc=[sp.factor(f) for f,m in facs if unitary_exact(f)[0]]
        cycpart.append((b,cyc))
    print("D_b, b=1..%d: 'monomial x cyclotomics' EXACTLY for b<=2 and never for b>=3 : %d mismatches"%(BM,nd))
    print("  cyclotomic part of D_b (exact factorisation):")
    for b,c in cycpart:
        if c: print("     b=%2d  (b mod 6 = %d) :"%(b,b%6), c)
    allcyc=all((not c) or (c==[t**2-t+1] or c==[sp.factor(t**2-t+1)]) for b,c in cycpart if b>=3)
    print("  every cyclotomic factor of D_b (b>=3) is Phi_6 :", allcyc)
    print("  b with a cyclotomic factor:", [b for b,c in cycpart if c], " = {b : b = 2 mod 6}?",
          [b for b,c in cycpart if c]==[b for b in range(1,BM+1) if b%6==2])
    # E_b
    ne=0
    for b in range(1,BM+1):
        ok,bad=unitary_exact(t**b-t**(b-1)+2)
        expect=(b==1)
        if ok!=expect: ne+=1; print("  MISMATCH E_%d ok=%s expect=%s"%(b,ok,expect))
    print("E_b, b=1..%d: 'monomial x cyclotomics' EXACTLY for b=1 and never for b>=2 : %d mismatches"%(BM,ne))
    # off-diagonal cells of Theorem C
    no=0; tot=0
    for b in range(1,31):
        cells=[(t-1)*t**(b-1)]+[(t-1)*t**(b-1)+(t-1)*t**(b-m-1) for m in range(1,b)]
        for f in cells:
            tot+=1; ok,bad=unitary_exact(f)
            if not ok: no+=1; print("  OFF-DIAG NOT UNITARY b=%d  %s  %s"%(b,sp.factor(f),bad))
    print("Theorem C off-diagonal cells (b<=30): %d cells, all 'monomial x cyclotomics' : %d failures"%(tot,no))
