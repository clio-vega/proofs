"""Case B1 elementary inequality g(j)=j-2s2(j)-4 v2(j(j-1)(j-2))+8 > 0 for j>=8, j not in {8,10}."""
import math
def v2(n):
    n=abs(n);k=0
    while n%2==0:n//=2;k+=1
    return k
def s2(n):return bin(n).count('1')
def g(j):return j-2*s2(j)-4*v2(j*(j-1)*(j-2))+8
bad=[j for j in range(8,4000) if g(j)<=0]
print("j>=8 with g(j)<=0 :",bad,"   (== {8,10})")
print("g(8..20):",[g(j) for j in range(8,21)])
thr=[j for j in range(32,100000) if j-6*math.log2(j)+2<=0]
print("j>=32 with j-6 log2 j+2<=0 :",thr,"(empty => analytic clears j>=32)")
