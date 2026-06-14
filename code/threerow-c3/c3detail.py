from math import comb
def Qf(av,bv,jv):
    H=( (av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4) )
    return (av-1)*(bv-2)*H-720*comb(jv,6)
def v2(n):
    n=int(round(n))
    if n==0:return 10**9
    c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)
def Delta(av,bv,jv):
    return jv-2*s2(jv)+2*(v2C(av+4,jv)+v2C(bv+3,jv)+v2(Qf(av,bv,jv))-v2(Qf(av,bv,0)))

# a EVEN: print Delta(1),Delta(2),Delta(4),Delta(6) by (a mod4,b mod4)
print("=== a EVEN: Delta(j) by (a%4,b%4) [b odd] ===")
import collections
tab=collections.defaultdict(lambda:collections.defaultdict(set))
for mv in range(4,60):
    n=2*mv
    for av in range(2,n+1,2):
        for bv in range(3,av+1,2):
            if n-av-bv!=3:continue
            for jv in [1,2,3,4,5,6]:
                if jv<=bv:
                    tab[(av%4,bv%4)][jv].add(Delta(av,bv,jv))
for cls in sorted(tab):
    row={jv:sorted(tab[cls][jv]) for jv in [1,2,3,4,5,6] if tab[cls][jv]}
    print(" (a%4,b%4)=",cls,row)

print("\n=== a ODD: forced descent Delta(1),Delta(2),Delta(3) and relative tildeDelta(j)=Delta(j)-Delta(3) ===")
tab2=collections.defaultdict(lambda:collections.defaultdict(set))
for mv in range(4,60):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            d3=Delta(av,bv,3)
            for jv in [1,2,3,5,7,9]:
                if jv<=bv:
                    val=Delta(av,bv,jv) if jv<=3 else Delta(av,bv,jv)-d3
                    tab2[(av%4,bv%4)][jv].add(val)
for cls in sorted(tab2):
    row={jv:sorted(tab2[cls][jv]) for jv in [1,2,3,5,7,9] if tab2[cls][jv]}
    print(" (a%4,b%4)=",cls,"  D1,D2,D3 then tilde(5,7,9):",row)
