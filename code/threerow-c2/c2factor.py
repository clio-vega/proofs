import sympy as sp
a,b,j,m,N = sp.symbols('a b j m N', integer=True)
# base binomial C(N, B) with B=b-j; N = a+b+2-2j  (since a+b=2m-2, N=2(m-j))
B = b-j
Nval = a+b+2-2*j

def ratio_same_N(d):
    # C(N, B+d)/C(N,B) as rational in a,b,j (N substituted)
    # C(N,k+1)/C(N,k) = (N-k)/(k+1)
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

def ratio_lowerN(drop, d):
    # C(N-drop, B+d)/C(N,B). First C(N-drop, B+d)/C(N,B+d) then * C(N,B+d)/C(N,B)
    # C(N-1,k)/C(N,k) = (N-k)/N ; iterate for drop
    r = ratio_same_N(d)
    k = B+d
    cur = sp.Integer(1)
    NN = N
    for s in range(drop):
        cur *= (NN-1-k)/(NN-1)   # C(NN-1,k)/C(NN,k) = (NN-k)/NN? wait
    # careful: C(M-1,k)/C(M,k) = (M-k)/M
    cur = sp.Integer(1); M = N
    for s in range(drop):
        cur *= (M-k)/M
        M = M-1
    return r*cur

# Now express each D_j term as (rational)*C(N,B). Need each term's binomials with their A=p-j.
# We'll build rho for each D_j call.
def rho_Dj0(p,q):
    A=p-j   # C(N, A) -> offset A-B
    d = sp.simplify(A-B)
    return ratio_same_N(int(d))
def rho_Dj1(p,q):
    A=p-j
    # N*C(N-1,A) + j*C(N,A) + j*C(N,A+1)
    d=int(sp.simplify(A-B))
    return N*ratio_lowerN(1,d) + j*ratio_same_N(d) + j*ratio_same_N(d+1)
def rho_Dj2(p,q):
    A=p-j
    d=int(sp.simplify(A-B))
    # t00 = C(N,2)*C(N-2,A) = N(N-1)/2 * C(N-2,A)
    t00 = N*(N-1)/2 * ratio_lowerN(2,d)
    # t10 = j*trinom(N;A,q-j+1,1) = j*(q-j+1)*C(N,A)   since trinom/.. ; verify: trinom(N,A,Y,1)=N!/(A!Y!1!) with A+Y+1=N => = C(N,A)*C(N-A,1)? = C(N,A)*(N-A). but N-A=Y+1=q-j+1. ok = (q-j+1)*C(N,A)
    t10 = j*(q-j+1)*ratio_same_N(d)
    # t01 = j*trinom(N;A+1,q-j,1)= j*(q-j)*C(N,A+1)? trinom(N,A+1,q-j,1), (A+1)+(q-j)+1=N => N-(A+1)=q-j. =C(N,A+1)*(q-j)
    t01 = j*(q-j)*ratio_same_N(d+1)
    # t20 = C(j,2)*trinom(N;A,q-j+2,0)=C(j,2)*C(N,A)
    t20 = sp.binomial(j,2)*ratio_same_N(d)
    # t02 = C(j,2)*C(N,A+2)
    t02 = sp.binomial(j,2)*ratio_same_N(d+2)
    # t11 = j(j-1)*C(N,A+1)
    t11 = j*(j-1)*ratio_same_N(d+1)
    return t00+t10+t01+t20+t02+t11

# M_j/C(N,B):
rho = ( rho_Dj2(a,b) - rho_Dj2(a+1,b-1) - rho_Dj1(a,b+1) - rho_Dj0(a+2,b)
        + rho_Dj0(a+1,b+1) + rho_Dj2(a+2,b-1) )
rho = sp.simplify(rho.subs(N, Nval))
print("M_j / C(N,b-j) =")
print(sp.factor(rho))
