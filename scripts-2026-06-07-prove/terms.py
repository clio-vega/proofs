import math
def v2(n):
    if n==0: return 10**9
    v=0
    while n%2==0: n//=2; v+=1
    return v
def C(n,k):
    if k<0 or k>n: return 0
    return math.comb(n,k)
def imv2(j):
    # v2 of Im((1+i)^j); Im=0 if j%4==0
    if j%4==0: return None
    return j//2
def term_v2(b,m,j):
    # term_j = C(m,j)*Im((1+i)^j)*(T(b-j)-T(b-1-j)), T(l)=[u^l](1+u^2)^{m-j}=C(m-j,l/2) if l even
    g=imv2(j)
    if g is None: return None
    # the binomial part: exactly one of b-j,b-1-j even
    l=b-j
    if l%2==0: beta=l//2
    else: beta=(b-1-j)//2
    binom=C(m-j,beta)
    if C(m,j)==0 or binom==0: return None
    return v2(C(m,j))+g+v2(binom)

for b in [5,8,9,12,13,16,17,20,21,24]:
    R=(b-1)//2
    bad=[]
    for m in range(b, b+80):
        t1=term_v2(b,m,1)
        worst=None
        for j in range(2, b+1):
            tj=term_v2(b,m,j)
            if tj is None: continue
            if worst is None or tj<worst: worst=tj
        if worst is not None and worst<=t1:
            bad.append((m,t1,worst))
    print(f"b={b:2d} b%4={b%4}: term1 dominates each j>=2 ? {'YES' if not bad else 'NO at '+str(bad[:5])}")
