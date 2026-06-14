from math import comb
def v2(n):
    n=int(n)
    if n==0: return 10**9
    r=0
    while n%2==0: n//=2; r+=1
    return r
def C(n,k): return comb(n,k) if 0<=k<=n else 0
def fall(n,k):
    p=1
    for i in range(k): p*=(n-i)
    return p
# Terms of C(b+3,j)*H(j) via subset identity. a even, any b.
# T_A = A*C(b+3,j); A=(a+3)(a+4)(b+2)(b+3)
# t2 = -6(a+3)(b+2)(b+3)C(b+2,j-1)
# t3 = -3(ab+a+2b)(b+3)(b+2)C(b+1,j-2)
# t4 = 6(b+3)(b+2)(b+1)C(b,j-3)
# t5 = 3(b+3)(b+2)(b+1)b C(b-1,j-4)
bad=0; tested=0; details={}
for a in range(2,80,2):
    for b in range(2,a+1):
        for j in range(0, b+4):
            tgt=v2(b+3)+1
            terms=[
              (a+3)*(a+4)*(b+2)*(b+3)*C(b+3,j),
              -6*(a+3)*(b+2)*(b+3)*C(b+2,j-1),
              -3*(a*b+a+2*b)*(b+3)*(b+2)*C(b+1,j-2),
              6*(b+3)*(b+2)*(b+1)*C(b,j-3),
              3*(b+3)*(b+2)*(b+1)*b*C(b-1,j-4),
            ]
            for ti,val in enumerate(terms):
                tested+=1
                if v2(val)<tgt:
                    bad+=1; details[ti]=details.get(ti,0)+1
print(f"term-by-term v2>=v2(b+3)+1 : failures {bad} of {tested}  per-term:{details}")
# also confirm sum of terms == C(b+3,j)*H(j)
def H(a,b,j):
    return (a+3)*(a+4)*(b+2)*(b+3) -6*(a+3)*(b+2)*C(j,1) -6*(a*b+a+2*b)*C(j,2) +36*C(j,3)+72*C(j,4)
bad2=0
for a in range(2,40,2):
    for b in range(2,a+1):
        for j in range(0,b+4):
            terms=[(a+3)*(a+4)*(b+2)*(b+3)*C(b+3,j),
              -6*(a+3)*(b+2)*(b+3)*C(b+2,j-1),
              -3*(a*b+a+2*b)*(b+3)*(b+2)*C(b+1,j-2),
              6*(b+3)*(b+2)*(b+1)*C(b,j-3),
              3*(b+3)*(b+2)*(b+1)*b*C(b-1,j-4)]
            if sum(terms)!=C(b+3,j)*H(a,b,j): bad2+=1
print("identity sum(terms)==C(b+3,j)H(j) failures:",bad2)
