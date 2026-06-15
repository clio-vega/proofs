def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def v2prod(lo,hi):
    s=0
    for t in range(lo,hi+1): s+=v2(t)
    return s
bad=0
for Q in range(4,80):
    for P in range(0,2000):
        if v2prod(P+1,P+Q)-1 < v2(P+Q-1): bad+=1
print("PRODUCT LEMMA (Q>=4): bad=",bad)
b3=0
for P in range(0,300):
    if v2prod(P+1,P+3)-1 < v2(P+3-1): b3+=1
print("Q=3 failures:",b3)
