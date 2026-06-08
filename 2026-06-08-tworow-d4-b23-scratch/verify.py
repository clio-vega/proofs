import sympy as sp
import random

def v2(n):
    n = int(n)
    if n == 0: return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

def C(n, k):
    if k < 0 or n < 0 or k > n: return 0
    return int(sp.binomial(n, k))

def Imii(j):
    return int(sp.im(sp.expand((1+sp.I)**j)))

def tau(j, b, m):
    def T(l):
        if l < 0 or l % 2 != 0: return 0
        return C(m-j, l//2)
    return int(C(m, j) * Imii(j) * (T(b-j) - T(b-1-j)))

def Ib(b, m):
    return sum(tau(j, b, m) for j in range(1, b+1))

# ---- Claim package ----
# R = (b-2)//2.  a := v2(tau_1).
# (P1) v2(tau_j) >= a for all j.
# Ties (j>=2) occur only for j in {2,3}, and exactly iff m odd.
# (D-forms) for j with Im((1+i)^j)!=0:  h=j//2
#    j odd : D0 = v2(C(R, h))
#    j even: D0 = v2(C(R, h-1)) - v2(h)
#    and  v2(tau_j) - a  ==  D0 + S_j,  S_j = sum_{t=1}^{h} v2(m-R-t)
# (V) v2(Ib) == a == v2((m)_{R+1}) - v2(R!).

def Spow(b, m, j):
    R = (b-2)//2
    h = j//2
    return sum(v2(m-R-t) for t in range(1, h+1))

def D0_formula(b, j):
    R = (b-2)//2
    h = j//2
    if j % 2 == 1:
        return v2(C(R, h))
    else:
        return v2(C(R, h-1)) - v2(h)

def v2_falling(m, k):  # v2 of (m)_k = m(m-1)...(m-k+1)
    return sum(v2(m-t) for t in range(0, k))

random.seed(1)
problems = 0
checks = 0
for b in range(4, 65, 4):       # b ≡ 0 mod 4
    R = (b-2)//2
    assert R % 2 == 1 and (R+1) % 2 == 0
    ms = list(range(b, b+20)) + [random.randint(b, 10**6) for _ in range(40)]
    for m in ms:
        a = v2(tau(1, b, m))
        # (V) RHS
        rhs = v2_falling(m, R+1) - v2(sp.factorial(R))
        assert a == rhs, ("a vs rhs", b, m, a, rhs)
        I = Ib(b, m)
        # (P1)+(V): check per term
        ties = []
        for j in range(1, b+1):
            t = tau(j, b, m)
            if t == 0:
                continue
            vj = v2(t)
            # (P1)
            if vj < a:
                print("DIP BELOW", b, m, j, vj, a); problems += 1
            # D-form check
            df = D0_formula(b, j) + Spow(b, m, j)
            if vj - a != df:
                print("DFORM MISMATCH", b, m, j, vj-a, df); problems += 1
            if vj == a:
                ties.append(j)
            checks += 1
        # tie structure
        expected_ties = [1,2,3] if m % 2 == 1 else [1]
        if ties != expected_ties:
            print("TIE STRUCT", b, m, "got", ties, "exp", expected_ties); problems += 1
        # (V)
        if v2(I) != a:
            print("V FAIL", b, m, v2(I), a); problems += 1
        # law
        if I == 0:
            print("LAW FAIL Ib=0", b, m); problems += 1

print(f"checks={checks} problems={problems}")
print("ALL OK" if problems == 0 else "PROBLEMS FOUND")
