from math import comb
import collections
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def Hf(av,bv,jv):
    return ((av+3)*(av+4)*(bv+2)*(bv+3)-6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4))
def Qf(av,bv,jv): return (av-1)*(bv-2)*Hf(av,bv,jv)-720*comb(jv,6)
def td(av,bv,jv):
    vP0=v2(Qf(av,bv,0))
    U=v2(comb(av+4,jv)*comb(bv+3,jv)*Qf(av,bv,jv))-vP0
    return jv+3-2*s2(jv)+2*U

for jv in (5,7,9):
    tab=collections.defaultdict(set)
    for mv in range(9,80):
        n=2*mv
        for av in range(1,n+1,2):           # a odd
            for bv in range(4,av+1,2):      # b even, b>=4
                if n-av-bv!=3:continue
                if bv<jv: continue          # interior j<=b
                tab[(av%4,bv%4)].add(td(av,bv,jv))
    print(f"j={jv}: td by (a%4,b%4):", {k:sorted(v) for k,v in sorted(tab.items())})

# also: confirm all-or-nothing - when a=1,b=2 mod4 and all of 5,7,9 interior (b>=9), all three td=0?
print("\nall-or-nothing check (b>=9):")
bad=0
for mv in range(9,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3 or bv<9:continue
            tds=[td(av,bv,j) for j in (5,7,9)]
            iscond=(av%4==1 and bv%4==2)
            allzero=all(t==0 for t in tds)
            anyzero=any(t==0 for t in tds)
            if allzero!=iscond or anyzero!=iscond: bad+=1
print("mismatch (all-zero<=>a1b2mod4, and any<=>all):",bad)
