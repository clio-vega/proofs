import sympy as sp
a_,b_=sp.symbols('a b',integer=True)
G2 = a_**2*b_**2+7*a_**2*b_+12*a_**2+9*a_*b_**2+31*a_*b_+28*a_+20*b_**2+28*b_+8
G2f=sp.lambdify((a_,b_),G2,'math')
def v2(n):
    n=abs(int(n))
    if n==0:return 99
    k=0
    while n%2==0:n//=2;k+=1
    return k
# T(2)=v2(G2)-2; tie (Delta(2)=0) iff v2(G2)=2.
from collections import defaultdict
tie_by_mod4=defaultdict(lambda:[0,0])  # (a%4,b%4)->[tie,total]
minv=99
for a in range(4,260):
    for b in range(4,a+1):
        if (a+b)%2: continue
        v=v2(G2f(a,b)); minv=min(minv,v)
        tie = (v==2)
        key=(a%4,b%4)
        tie_by_mod4[key][1]+=1
        if tie: tie_by_mod4[key][0]+=1
print("min v2(G2)=",minv,"(need >=2)")
print("tie (v2(G2)=2) by (a%4,b%4):")
for k in sorted(tie_by_mod4):
    t,tot=tie_by_mod4[k]
    print(f"  {k}: {t}/{tot}  {'ALWAYS TIE' if t==tot else ('NEVER' if t==0 else 'MIXED')}")
