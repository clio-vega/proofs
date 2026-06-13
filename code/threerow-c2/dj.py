from mn import Mj, fdim
from math import comb, factorial
from itertools import permutations

def trinom(N,a,b,c):
    if a<0 or b<0 or c<0 or a+b+c!=N: return 0
    return factorial(N)//(factorial(a)*factorial(b)*factorial(c))

def Dj(p,q,r,j,m):
    # [s^p t^q u^r] e1^{2(m-j)} e2^j
    if p<0 or q<0 or r<0: return 0
    if p+q+r != 2*m: return 0
    N=2*(m-j)
    tot=0
    for i in range(j+1):
        for k in range(j-i+1):
            l=j-i-k
            # s exp: i+k, t exp: i+l, u exp: k+l, plus e1 part
            coef = factorial(j)//(factorial(i)*factorial(k)*factorial(l))
            tot += coef * trinom(N, p-(i+k), q-(i+l), r-(k+l))
    return tot

def Mj_jt(lam,j,m):
    a,b,c=lam
    L=[a-1,b-2,c-3]
    tot=0
    for sg,sigma in [(1,(1,2,3)),(-1,(1,3,2)),(-1,(2,1,3)),(1,(2,3,1)),(1,(3,1,2)),(-1,(3,2,1))]:
        idx=[L[0]+sigma[0], L[1]+sigma[1], L[2]+sigma[2]]
        tot += sg*Dj(idx[0],idx[1],idx[2],j,m)
    return tot

bad=0
for m in range(2,9):
    n=2*m
    for a in range(1,n+1):
        for b in range(1,a+1):
            c=n-a-b
            if c<1 or c>b: continue
            for j in range(0,m+1):
                if Mj((a,b,c),j,m)!=Mj_jt((a,b,c),j,m):
                    bad+=1
                    if bad<10: print("MISMATCH",(a,b,c),m,j,Mj((a,b,c),j,m),Mj_jt((a,b,c),j,m))
print("D_j/JT formula check, mismatches =",bad)
