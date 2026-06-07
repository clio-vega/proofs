import math
def v2(n):
    n=int(n)
    if n==0: return 10**9
    v=0
    while n%2==0: n//=2; v+=1
    return v
def vfact(n): return v2(math.factorial(n))

def check_star(b):
    R=(b-1)//2
    bad=[]
    for j in range(2,b+1):
        if j%4==0: continue            # Im=0, term absent
        dj=(j-1)//2
        sj=R+1+dj
        if sj<j: continue              # C(sj,j)=0, term absent
        lhs=v2(math.comb(sj,j))+j//2
        rhs=1+(vfact(sj)-vfact(R+1))
        if lhs<rhs: bad.append((j,lhs,rhs))
    return bad

print("b=1 mod4: verify (*) for all valid j")
for b in range(5,120,4):  # b=5,9,13,...
    bad=check_star(b)
    print(f"  b={b:3d} R={(b-1)//2:3d}: {'(*) holds' if not bad else 'FAILS '+str(bad[:6])}")
