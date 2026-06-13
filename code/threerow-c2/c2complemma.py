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

# Test sufficient condition: v2 Q(j) >= v2(a(b-1)) + v2(j) + 1   (for 1<=j<=b)
suff_bad=0; suff_ex=[]
# Also test the REAL compensation lemma T(j)>=1-v2(j)
lem_bad=0
# And B(j)>=1 off {2,4}
B_bad=0
for A in range(3,300):
    for B in range(2,A+1):
        if (A+B)%2: continue
        for j in range(1,B+1):
            g=v2(A*(B-1))
            if v2(Q(A,B,j)) < g+v2(j)+1:
                suff_bad+=1
                if len(suff_ex)<8: suff_ex.append((A,B,j,v2(Q(A,B,j)),g+v2(j)+1))
            T=vC(A+3,j)+vC(B+2,j)+v2(Q(A,B,j))-v2(Q(A,B,0))
            if T < 1-v2(j): lem_bad+=1
for j in range(1,2000):
    if j in (2,4): continue
    if j-2*s2(j)+2-2*v2(j) < 1: B_bad+=1
print("sufficient cond v2Q(j)>=v2(a(b-1))+v2(j)+1 FAILURES:",suff_bad, suff_ex[:8])
print("REAL Comp Lemma T(j)>=1-v2(j) FAILURES:",lem_bad)
print("B(j)>=1 off {2,4} FAILURES (j<2000):",B_bad)
