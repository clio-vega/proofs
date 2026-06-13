import sys
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def vC(n,k):
    if k<0 or k>n:return 10**9
    return v2(comb(n,k))
def Q(a,b,j): return a*(b-1)*((a+3)*(b+2)-2*j*j)+j*(j-1)*(j-2)*(j-3)
def Delta(a,b,j): return j-2*s2(j)+2*vC(a+3,j)+2*vC(b+2,j)+2*(v2(Q(a,b,j))-v2(Q(a,b,0)))

neg2=neg4=0
from collections import defaultdict
z2=defaultdict(set); z4=defaultdict(set); both0=0
for A in range(4,400):
    for B in range(4,A+1):
        if (A+B)%2: continue
        d2=Delta(A,B,2); d4=Delta(A,B,4)
        if d2<0: neg2+=1
        if d4<0: neg4+=1
        z2[(A%4,B%4)].add(d2==0)
        z4[(A%4,B%4)].add(d4==0)
        if d2==0 and d4==0: both0+=1
print("Delta(2)<0 count:",neg2," Delta(4)<0 count:",neg4)
print("both Delta(2)=0 and Delta(4)=0:",both0)
print("\n{Delta(2)==0} by (a%4,b%4):")
for k in sorted(z2): print("  ",k,"-> zero?",z2[k])
print("{Delta(4)==0} by (a%4,b%4):")
for k in sorted(z4): print("  ",k,"-> zero?",z4[k])
