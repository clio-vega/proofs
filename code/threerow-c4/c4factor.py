import sympy as sp
a,b,j,m,N = sp.symbols('a b j m N', integer=True)

# c=4:  a+b+4 = 2m  => a+b = 2m-4,  N = 2(m-j) = a+b+4-2j
B = b-j
Nval = a+b+4-2*j

def ratio_same_N(d):
    # C(N, B+d)/C(N,B) as rational, N symbolic.  d a Python int.
    expr = sp.Integer(1)
    if d>=0:
        k=B
        for _ in range(d):
            expr*= (N-k)/(k+1); k=k+1
    else:
        k=B-1
        for _ in range(-d):
            expr*= (k+1)/(N-k); k=k-1
    return expr

def falling_binom(top_expr, r):
    if r<0: return sp.Integer(0)
    if r==0: return sp.Integer(1)
    num=sp.Integer(1)
    for t in range(r):
        num*= (top_expr - t)
    return num/sp.factorial(r)

def rho_Dj(q,r,qb):
    # D_j(p,q,r)/C(N,b-j) referenced to t-coordinate. qb = q-b.
    # D_j = sum_w C(j,w) sum_l C(w,l) C(N, q-j+w-l) C(N-(q-j+w-l), r-w)
    total=sp.Integer(0)
    for w in range(r+1):
        for l in range(w+1):
            coef = sp.binomial(j,w)*sp.binomial(w,l)
            off = qb + (w-l)
            rr = ratio_same_N(off)
            te1 = q - j + (w-l)
            fb = falling_binom(N - te1, r-w)
            total += coef*rr*fb
    return total

# c=4 expansion:
# M_j = D(a,b,4) - D(a+1,b-1,4) - D(a,b+1,3) - D(a+2,b,2) + D(a+1,b+1,2) + D(a+2,b-1,3)
rho = ( rho_Dj(b,   4,  0)
      - rho_Dj(b-1, 4, -1)
      - rho_Dj(b+1, 3,  1)
      - rho_Dj(b,   2,  0)
      + rho_Dj(b+1, 2,  1)
      + rho_Dj(b-1, 3, -1) )
rho = sp.simplify(rho.subs(N, Nval))
print("M_j / C(N,b-j) =")
num,den = sp.fraction(sp.together(rho))
num=sp.expand(num); den=sp.expand(den)
print("NUM=",sp.factor(num))
print("DEN=",sp.factor(den))
