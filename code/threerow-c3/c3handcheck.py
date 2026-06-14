from math import comb
def Hf(av,bv,jv):
    return ( (av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4) )
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
def Delta(av,bv,jv):
    return jv-2*s2(jv)+2*(v2C(av+4,jv)+v2C(bv+3,jv)+v2(Qf(av,bv,jv))-v2(Qf(av,bv,0)))

c1=c2=c3=c4=c5=0
for mv in range(4,60):
    n=2*mv
    for av in range(2,n+1,2):           # a even
        for bv in range(3,av+1,2):      # b odd
            if n-av-bv!=3:continue
            # claim: Delta(1)=1
            if Delta(av,bv,1)!=1: c1+=1
            # claim: Delta(2)=2*(v2H(2)-2)
            if bv>=2 and Delta(av,bv,2)!=2*(v2(Hf(av,bv,2))-2): c2+=1
            # claim: Delta(2)=0 <=> a%4==2 and b%4==3
            if bv>=2:
                tie2=(Delta(av,bv,2)==0)
                cong2=(av%4==2 and bv%4==3)
                if tie2!=cong2: c3+=1
            # claim: Delta(4)=0 <=> a%4==0 and b%4==1
            if bv>=4:
                tie4=(Delta(av,bv,4)==0)
                cong4=(av%4==0 and bv%4==1)
                if tie4!=cong4: c4+=1
            # claim v2H(2)=2 <=> a%4==2,b%4==3  (the Delta(2) tie condition)
            if bv>=2:
                if (v2(Hf(av,bv,2))==2)!=(av%4==2 and bv%4==3): c5+=1
print("a even hand-claim failures:")
print("  Delta(1)=1:",c1)
print("  Delta(2)=2(v2H(2)-2):",c2)
print("  Delta(2)=0 <=> a=2,b=3 mod4:",c3)
print("  Delta(4)=0 <=> a=0,b=1 mod4:",c4)
print("  v2H(2)=2 <=> a=2,b=3 mod4:",c5)
