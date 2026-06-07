import math, random
def v2(n):
    n=int(n)
    if n==0: return 10**9
    v=0
    while n%2==0: n//=2; v+=1
    return v
def C(n,k):
    if k<0 or k>n: return 0
    return math.comb(n,k)
def imv2(j):
    if j%4==0: return None
    return j//2
def term_v2(b,mm,j):
    g=imv2(j)
    if g is None: return None
    l=b-j
    beta = l//2 if l%2==0 else (b-1-j)//2
    cm=C(mm,j); bn=C(mm-j,beta)
    if cm==0 or bn==0: return None
    return v2(cm)+g+v2(bn)

ok=True
for b in [5,9,13,17,21,25,29,33,37,41]:   # b = 1 mod 4
    R=(b-1)//2
    for _ in range(300):
        mm=random.randint(b, 100000)
        t1=term_v2(b,mm,1)
        for j in range(2,b+1):
            tj=term_v2(b,mm,j)
            if tj is None: continue
            h=j//2; dj=(j-1)//2
            Dj=v2(h+1)+v2(C(R+1,h+1))
            S=sum(v2(mm-R-t) for t in range(1,dj+1))
            pred = t1 + Dj + S
            if pred != tj:
                ok=False
                print(f"MISMATCH b={b} m={mm} j={j}: actual {tj} pred {pred} (Dj={Dj} S={S})")
            if tj <= t1:
                ok=False
                print(f"NONDOM b={b} m={mm} j={j}: tj={tj} t1={t1}")
print("ALL DECOMPOSITION & STRICT DOMINATION OK" if ok else "PROBLEMS FOUND")
