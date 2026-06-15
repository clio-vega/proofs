from math import comb, factorial

def v2(n):
    if n == 0: return float('inf')
    k = 0
    while n % 2 == 0:
        n//=2; k+=1
    return k

def descfact(j, k):
    p = 1
    for i in range(k):
        p *= (j - i)
    return p

# (1) subset identity: C(F+2c-1, j)*C(j,2c) = C(F+2c-1,2c)*C(F-1, j-2c)
# valid for 2c <= j <= F+2c-1, F>=1 (so F-1>=0)
bad1 = 0
for c in range(1, 7):
    for F in range(2, 202, 2):
        N = F + 2*c - 1
        for j in range(2*c, N+1):
            lhs = comb(N, j) * comb(j, 2*c)
            rhs = comb(N, 2*c) * comb(F-1, j-2*c)
            if lhs != rhs:
                bad1 += 1
print("subset identity mismatches:", bad1)

# (2) central-tip: j^(2c) = (2c)! * C(j,2c)  for all j
bad2 = 0
for c in range(1, 7):
    for j in range(0, 200):
        if descfact(j, 2*c) != factorial(2*c) * comb(j, 2*c):
            bad2 += 1
print("central-tip mismatches:", bad2)

# (3) c=3 anchor identity: v2(C(F+5,6)) + v2(6!) = sum_{i=0}^{5} v2(F+i)
bad3 = 0
for F in range(2, 80, 2):
    lhs = v2(comb(F+5, 6)) + v2(factorial(6))
    rhs = sum(v2(F+i) for i in range(6))
    if lhs != rhs:
        bad3 += 1
print("c=3 anchor mismatches:", bad3)
