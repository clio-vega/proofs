from sympy import factorial
def v2(n):
    n=int(n)
    if n==0: return 10**9
    r=0
    while n%2==0: n//=2; r+=1
    return r
from math import comb
def s2(n):
    return bin(n).count('1')
def beta(c):
    return (c-1)+v2(factorial(c-1))
def falling(j,k):
    p=1
    for i in range(k): p*=(j-i)
    return p

# 1) beta formula vs closed form 2(c-1)-s2(c-1)
print("c, beta=(c-1)+v2((c-1)!), 2(c-1)-s2(c-1)")
for c in range(1,9):
    print(c, beta(c), 2*(c-1)-s2(c-1))

# 2) Verify NL_c directly and tightness
print("\nNL_c check (min slack should be 0 => beta tight):")
for c in range(1,6):
    b=beta(c)
    worst=10**9; bad=0
    for F in range(2, 200, 2):
        for j in range(2*c, F+2*c):
            lhs=v2(comb(F+2*c-1,j))+v2(falling(j,2*c))
            rhs=v2(F)+b
            slack=lhs-rhs
            if slack<0: bad+=1
            worst=min(worst,slack)
    print(f"c={c} beta={b} min_slack={worst} failures={bad}")
