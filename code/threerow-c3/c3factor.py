import sympy as sp
a,b,j,m,N = sp.symbols('a b j m N', integer=True)

# For c=3: a+b+3 = 2m  => a+b = 2m-3, N = 2(m-j) = a+b+3-2j
B = b-j
Nval = a+b+3-2*j

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

def rho_Dj(p,q,r,qb):
    # D_j(p,q,r)/C(N,b-j). qb = q-b (small int, in {-1,0,1}).
    # D_j = sum_w C(j,w) sum_l C(w,l) C(N, q-j+w-l) C(N-(q-j+w-l), r-w)
    total=sp.Integer(0)
    for w in range(r+1):
        for l in range(w+1):
            coef = sp.binomial(j,w)*sp.binomial(w,l)
            off = qb + (w-l)            # t_e1 - B
            rr = ratio_same_N(off)
            te1 = q - j + (w-l)         # symbolic ok via Nval later
            fb = falling_binom(N - te1, r-w)
            total += coef*rr*fb
    return total

# c=3 expansion: M_j = D(a,b,3) - D(a,b+1,2) - D(a+1,b-1,3) + D(a+1,b+1,1) + D(a+2,b-1,2) - D(a+2,b,1)
# qb = q - b for each:
rho = ( rho_Dj(a,   b,   3,  0)
      - rho_Dj(a,   b+1, 2,  1)
      - rho_Dj(a+1, b-1, 3, -1)
      + rho_Dj(a+1, b+1, 1,  1)
      + rho_Dj(a+2, b-1, 2, -1)
      - rho_Dj(a+2, b,   1,  0) )
rho = sp.simplify(rho.subs(N, Nval))
print("M_j / C(N,b-j) =")
num,den = sp.fraction(sp.together(rho))
print("NUM=",sp.factor(num))
print("DEN=",sp.factor(den))
