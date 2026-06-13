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

# verify Delta(2)=2[v2(a+2)+v2(b+1)+v2(W-8)-2]
bad2=0
for A in range(4,400):
    for B in range(4,A+1):
        if (A+B)%2: continue
        W=(A+3)*(B+2)
        f=2*(v2(A+2)+v2(B+1)+v2(W-8)-2)
        if f!=Delta(A,B,2): bad2+=1
print("Delta(2) closed-form mismatches:",bad2)

# T(j) lower bound study
def T(a,b,j): return vC(a+3,j)+vC(b+2,j)+v2(Q(a,b,j))-v2(Q(a,b,0))
minT={}  # by j
for A in range(4,260):
    for B in range(4,A+1):
        if (A+B)%2: continue
        for j in range(1,B+1):
            t=T(A,B,j)
            if j not in minT or t<minT[j][0]:
                minT[j]=(t,A,B)
print("\nmin T(j) by j (j: minT, witness a,b):")
for j in sorted(minT):
    print(f"  j={j}: minT={minT[j][0]}  at a,b={minT[j][1],minT[j][2]}")
