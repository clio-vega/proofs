"""Gap 2 (a-odd) exploration. Understand tildeDelta structure before proving.
tildeDelta(j) = j+3-2s2(j)+2U(j), U=v2C(a+4,j)+v2C(b+3,j)+v2Q3(j)-v2Q3(0).
a odd, b even, a>=b>=4, a+b=2m-3.  Interior 4<=j<=b. Offset j0=3.
"""
from math import comb
import collections

def Hf(av,bv,jv):
    return ((av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4))
def Qf(av,bv,jv):
    return (av-1)*(bv-2)*Hf(av,bv,jv)-720*comb(jv,6)
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)

def U(av,bv,jv):
    return v2C(av+4,jv)+v2C(bv+3,jv)+v2(Qf(av,bv,jv))-v2(Qf(av,bv,0))
def td(av,bv,jv):
    return jv+3-2*s2(jv)+2*U(av,bv,jv)

# 1. global: min tildeDelta and where 0 occurs
mins=10**9; fails=0
tie_js=collections.Counter()
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            for jv in range(4,bv+1):
                t=td(av,bv,jv)
                mins=min(mins,t)
                if t<0:fails+=1
                if t==0:tie_js[jv]+=1
print("min tildeDelta:",mins,"fails:",fails)
print("tie j-values:",dict(tie_js))

# 2. For tie shapes, show v2Q3(j)-v2Q3(0), v2C(a+4,j), v2C(b+3,j), s2 at j=5,7,9
print("\n-- tie shapes a=1,b=2 mod4, detail at j=5,7,9 --")
cnt=0
for mv in range(7,30):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            if av%4==1 and bv%4==2 and bv>=9:
                cnt+=1
                if cnt>6:continue
                row=[]
                for jv in [5,7,9]:
                    row.append((jv, v2C(av+4,jv),v2C(bv+3,jv),
                               v2(Qf(av,bv,jv))-v2(Qf(av,bv,0)), td(av,bv,jv)))
                print(f"a={av}(%4={av%4}) b={bv}(%4={bv%4}): "+
                      " | ".join(f"j={r[0]} vCa={r[1]} vCb={r[2]} dQ={r[3]} td={r[4]}" for r in row))
