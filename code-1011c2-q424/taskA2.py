import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1010-q410')
exec(open('taskA.py').read().split('print("="*78)')[0])

print("="*78); print("TASK A.2  SWEEP  n=2, alpha=(1,1), k=1..4, ell=0..3  (Rick's 16 pairs)"); print("="*78)
rows=[]
for k in range(1,5):
    for ell in range(0,4):
        h,w = 2*k, 2*ell+2
        nd  = nu_data(k,ell,(1,1))
        minnu = min(x[1] for x in nd) if nd else None
        c1  = c_alpha_M1(k,ell,(1,1))
        c2,A,B,TOT = c_alpha_M2(k,ell,(1,1))
        dom = (h>2) and (2<=w)
        ks  = sp.expand(sum(sign_hkko(la,k,w)*psi for la,nu,psi in nd))
        rows.append((k,ell,h,w,len(nd),minnu,sp.expand(c1),ordt1(c1),c2,dom,sp.expand(c1).subs(t,1)==c2, sp.simplify(ks-c1)==0))
print(f"{'k':>2}{'ell':>4}{'h':>3}{'w':>3}{'#T':>4}{'minnu':>6}  {'c_alpha(t) [M1]':<16}{'ord':>4}{'c(1) [M2]':>10}{'inDom':>7}{'M1=M2':>7}{'KorffOK':>8}")
for r in rows:
    print(f"{r[0]:>2}{r[1]:>4}{r[2]:>3}{r[3]:>3}{r[4]:>4}{r[5]:>6}  {str(r[6]):<16}{str(r[7]):>4}{r[8]:>10}{str(r[9]):>7}{str(r[10]):>7}{str(r[11]):>8}")
agree = sum(1 for r in rows if r[10]); print(f"  M1=M2 at t=1: {agree}/{len(rows)}")
print(f"  c_alpha in {{1-t, 1}}: {sum(1 for r in rows if r[6] in (1, sp.expand(1-t)))}/{len(rows)}  (Rick's reported set)")
print(f"  obstruction fires (ord < minnu) in-domain rows: "
      f"{sum(1 for r in rows if r[9] and isinstance(r[7],int) and r[7]<r[5])}"
      f"/{sum(1 for r in rows if r[9])}")
print(f"  Korff identity FAILS: {sum(1 for r in rows if not r[11])}/{len(rows)}")

print()
print("="*78); print("TASK A.3  CONTROL (i)  constant step => nu=0 => obstruction must VANISH"); print("="*78)
for (k,ell,alpha,why) in [(1,0,(2,),"h_cyc=2, alpha_1=2=h => theta=(1,1) CONSTANT"),
                          (1,1,(2,),"h_cyc=2, alpha_1=2=h => theta=(1,1) CONSTANT"),
                          (2,1,(4,),"h_cyc=4, alpha_1=4=h => theta=(1,1,1,1) CONSTANT"),
                          (2,1,(2,),"h_cyc=4, alpha_1=2 NOT in {0,h} => NON-constant")]:
    nd = nu_data(k,ell,alpha); h,w=2*k,2*ell+2
    minnu = min(x[1] for x in nd) if nd else None
    c1 = c_alpha_M1(k,ell,alpha); c2,A,B,TOT = c_alpha_M2(k,ell,alpha)
    ks = sp.expand(sum(sign_hkko(la,k,w)*psi for la,nu,psi in nd))
    print(f"  k={k} ell={ell} alpha={alpha} (h={h},w={w}): #T={len(nd)} minnu={minnu}  "
          f"c_alpha={sp.expand(c1)} ord={ordt1(c1)}  M2 c(1)={c2}  "
          f"obstruction vacuous(minnu=0)={minnu==0}  KorffIdentity={sp.simplify(ks-c1)==0}")
    print(f"      [{why}]")

print()
print("="*78); print("TASK A.4  CONTROL (ii)  ell=0 must give 1-t (ord 1), ell>=1 must give 1 (ord 0)"); print("="*78)
for k in range(1,5):
    r0 = [r for r in rows if r[0]==k and r[1]==0][0]
    r1 = [r for r in rows if r[0]==k and r[1]==1][0]
    print(f"  k={k}: ell=0 -> c={r0[6]} ord={r0[7]}   |   ell=1 -> c={r1[6]} ord={r1[7]}   "
          f"DISTINCT ORDERS = {r0[7]!=r1[7]}")

print()
print("="*78); print("TASK A.5  CONTROL (iii) REGRESSION: out-of-domain k=1, alpha=(1^n), n=2M"); print("="*78)
for ell in range(0,3):
    k=1; M=k+ell+1; n=2*M
    c2,A,B,TOT = c_alpha_M2(k,ell,tuple([1]*n))
    CM = sp.binomial(2*M,M)-sp.binomial(2*M,M-1)
    pred = sp.expand(CM-(2*M-1)*t**(M-1)+t**(M+1))
    print(f"  k=1 ell={ell} M={M} n={n}: M2 c(1)={c2}; Q405 closed form {pred} has c(1)={pred.subs(t,1)}  "
          f"MATCH={c2==pred.subs(t,1)}  (n={n} > w={2*ell+2} => OUT of Korff's domain)")
