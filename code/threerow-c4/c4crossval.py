import sympy as sp
a_,b_,j_=sp.symbols('a b j',integer=True)
H=(a_**3*b_**3+9*a_**3*b_**2+26*a_**3*b_+24*a_**3+12*a_**2*b_**3-4*a_**2*b_**2*j_**2-8*a_**2*b_**2*j_+108*a_**2*b_**2-12*a_**2*b_*j_**2-48*a_**2*b_*j_+312*a_**2*b_-8*a_**2*j_**2-64*a_**2*j_+288*a_**2+47*a_*b_**3-20*a_*b_**2*j_**2-64*a_*b_**2*j_+423*a_*b_**2+6*a_*b_*j_**4-12*a_*b_*j_**3-30*a_*b_*j_**2-384*a_*b_*j_+1222*a_*b_+24*a_*j_**3-40*a_*j_**2-488*a_*j_+1128*a_+60*b_**3-24*b_**2*j_**2-120*b_**2*j_+540*b_**2+6*b_*j_**4+12*b_*j_**3-42*b_*j_**2-696*b_*j_+1560*b_-4*j_**6+48*j_**5-244*j_**4+648*j_**3-664*j_**2-648*j_+1440)
P8=sp.prod([j_-t for t in range(8)])
Q4=sp.expand((a_-2)*(b_-3)*H+P8)
Hf=sp.lambdify((a_,b_,j_),H,'math'); Q4f=sp.lambdify((a_,b_,j_),Q4,'math')
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n):return bin(n).count('1')
def vfact(r):return sum(v2(t) for t in range(1,r+1))
def vfall(top,r):return 0 if r<=0 else sum(v2(top-t) for t in range(r))
def Delta(av,bv,jv):  # full Prop-2 Delta
    return jv-2*s2(jv)+2*(vfall(av+2,jv-3)+vfall(bv+1,jv-3)-2*vfact(jv)+v2(Q4f(av,bv,jv))-v2(av-2)-v2(bv-3))
def g(j):return j-2*s2(j)-4*v2(j*(j-1)*(j-2))+8

# For all shapes (a,b,4), m<=80, every interior j>=8: classify and check proven bound
bad=0;chk=0
caseA=caseB=0
for av in range(8,160):
    for bv in range(8,av+1):
        if (av+bv)%2:continue
        for jv in range(8,bv+1):
            chk+=1
            D=Delta(av,bv,jv)
            vheavy=v2(av-2)+v2(bv-3)+v2(Hf(av,bv,jv))
            vP8=vfall(jv,8)  # = v2 P8(j)
            if vheavy>=vP8:    # Case A
                caseA+=1
                ok = (D>=3)            # proven: Delta>=3
            else:              # Case B
                caseB+=1
                if jv in (8,10):
                    ok = (D>=2)        # finite check gives Delta>=2
                else:
                    ok = (D>=g(jv) and g(jv)>0)  # crude bound
            if D<=0 or not ok:
                bad+=1
                if bad<8:print("  PROBLEM",(av,bv,jv),"D=",D,"caseA" if vheavy>=vP8 else "caseB","g=",g(jv))
print(f"checked {chk} interior (j>=8) triples; CaseA={caseA} CaseB={caseB}; problems={bad}")
print("All Delta>0 and each within its proven bound." if bad==0 else "SOMETHING WRONG")

# Also: is N4 floor (v2H>=4) ever the *binding* constraint that the crude bound needs? show min over caseB nonspecial
mn=10**9
for av in range(8,160):
    for bv in range(8,av+1):
        if (av+bv)%2:continue
        for jv in range(8,bv+1):
            vheavy=v2(av-2)+v2(bv-3)+v2(Hf(av,bv,jv));vP8=vfall(jv,8)
            if vheavy<vP8 and jv not in(8,10):
                mn=min(mn,Delta(av,bv,jv)-g(jv))
print("min (Delta - g(j)) over CaseB nonspecial (>=0 confirms crude bound valid):",mn)
