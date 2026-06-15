"""Probe the TIP inequality. After absorbing C(j,6) into b-side:
  720 C(b+3,j) C(j,6) = (b+3)^(6) C(b-3,j-6),  v2((b+3)^(6)) = b2+v2(b)+b2p   (b2=v2(b-2),b2p=v2(b+2))
(star2''):  v2C(a+4,j) + v2(b) + v2C(b-3,j-6) >= a1+a3p + g(j),
  a1=v2(a-1), a3p=v2(a+3), g(j)=s2(j)-(j+3)/2.   For 6<=j<=b, a odd, b even.
Test star2'' and simpler variants to find the cleanest TRUE sub-lemma.
"""
from math import comb
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)

# candidate bounds, count failures
cands={
 "star2'': vCa + v2b + vC(b-3,j-6) >= a1+a3p+g":0,
 "A: vCa + v2b >= a1+a3p+g (drop vC(b-3))":0,
 "B: vCa >= a1+a3p+g-1 ":0,
 "C: vC(a+4,j)+vC(a-?,) ... use NL form vCa+vC(j,6) >= a1-1":0,
}
# also track margins
import collections
marg=collections.defaultdict(lambda:10**9)
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            a1=v2(av-1); a3p=v2(av+3)
            for jv in range(6,bv+1):
                g=s2(jv)-(jv+3)/2  # half-int allowed; compare 2x
                vCa=v2C(av+4,jv)
                # star2'': use 2x to stay integer
                L=2*(vCa+v2(bv)+v2C(bv-3,jv-6)); R=2*(a1+a3p)+2*s2(jv)-(jv+3)
                if L<R: cands["star2'': vCa + v2b + vC(b-3,j-6) >= a1+a3p+g"]+=1
                L2=2*(vCa+v2(bv));
                if L2<R: cands["A: vCa + v2b >= a1+a3p+g (drop vC(b-3))"]+=1
                L3=2*vCa+2;
                if L3<R: cands["B: vCa >= a1+a3p+g-1 "]+=1
                # NL form: vC(a+4,j)+vC(j,6) >= a1-1
                if vCa+v2C(jv,6) < a1-1: cands["C: vC(a+4,j)+vC(a-?,) ... use NL form vCa+vC(j,6) >= a1-1"]+=1
                marg[v2(jv)]=min(marg[v2(jv)], L-R)  # star2'' margin (x2)
for k,v in cands.items(): print(f"{v:5d} fails  | {k}")
print("star2'' margin (x2) by v2(j):",dict(sorted(marg.items())))
