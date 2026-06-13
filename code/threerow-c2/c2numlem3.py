def s2(n): return bin(n).count('1')
def v2(n):
    if n==0: return 10**9
    n=abs(n); k=0
    while n%2==0: n//=2;k+=1
    return k
def vC(n,k):
    if k<0 or k>n: return 10**9
    return s2(k)+s2(n-k)-s2(n)
def vP4(j):
    if j<4: return 10**9
    return v2(j)+v2(j-1)+v2(j-2)+v2(j-3)
bad=0; ex=[]
for F in range(2,3000,2):
    vF=v2(F); rhs=vF+1
    for j in range(4,F+4):
        if vC(F+3,j)+vP4(j) < rhs:
            bad+=1
            if len(ex)<15: ex.append((F,j,vF,vC(F+3,j),vP4(j)))
print("Number Lemma FAILURES (F<3000):",bad)
print("fail examples:",ex[:15])
