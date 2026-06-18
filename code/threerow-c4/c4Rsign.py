import sympy as sp
a_,b_,j_ = sp.symbols('a b j', integer=True)
H = (a_**3*b_**3 + 9*a_**3*b_**2 + 26*a_**3*b_ + 24*a_**3 + 12*a_**2*b_**3
 - 4*a_**2*b_**2*j_**2 - 8*a_**2*b_**2*j_ + 108*a_**2*b_**2 - 12*a_**2*b_*j_**2
 - 48*a_**2*b_*j_ + 312*a_**2*b_ - 8*a_**2*j_**2 - 64*a_**2*j_ + 288*a_**2
 + 47*a_*b_**3 - 20*a_*b_**2*j_**2 - 64*a_*b_**2*j_ + 423*a_*b_**2 + 6*a_*b_*j_**4
 - 12*a_*b_*j_**3 - 30*a_*b_*j_**2 - 384*a_*b_*j_ + 1222*a_*b_ + 24*a_*j_**3
 - 40*a_*j_**2 - 488*a_*j_ + 1128*a_ + 60*b_**3 - 24*b_**2*j_**2 - 120*b_**2*j_
 + 540*b_**2 + 6*b_*j_**4 + 12*b_*j_**3 - 42*b_*j_**2 - 696*b_*j_ + 1560*b_
 - 4*j_**6 + 48*j_**5 - 244*j_**4 + 648*j_**3 - 664*j_**2 - 648*j_ + 1440)
P8 = sp.prod([j_-t for t in range(8)])
Q4 = sp.expand((a_-2)*(b_-3)*H + P8)
Q4f = sp.lambdify((a_,b_,j_), Q4, 'math')

def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2; k+=1
    return k

# Scan R(j) = v2 Q4(j) - v2 Q4(0), find min over shapes per j.
minR={}
argmin={}
for av in range(4,200):
    for bv in range(4,av+1):
        if (av+bv)%2!=0: continue
        q0=Q4f(av,bv,0)
        v0=v2(q0)
        for jv in range(1,min(bv,40)+1):
            R=v2(Q4f(av,bv,jv))-v0
            if jv not in minR or R<minR[jv]:
                minR[jv]=R; argmin[jv]=(av,bv)
print("min R(j) per j (a,b<200 same parity):")
for jv in sorted(minR):
    print(f"  j={jv}: minR={minR[jv]}  at {argmin[jv]}")
