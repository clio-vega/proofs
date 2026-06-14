import sympy as sp
a,b,j=sp.symbols('a b j',integer=True)
# Q3 via structural identity, double-checked numerically
H=( (a+3)*(a+4)*(b+2)*(b+3) -6*(a+3)*(b+2)*sp.binomial(j,1)
    -6*(a*b+a+2*b)*sp.binomial(j,2)+36*sp.binomial(j,3)+72*sp.binomial(j,4) )
Q3id=(a-1)*(b-2)*H-720*sp.binomial(j,6)

# original Q3
K3=(a-1)*(a+3)*(a+4)*(b-2)*(b+2)*(b+3)
pert=(-j**6+15*j**5 + j**4*(3*a*b-6*a-3*b-79) + j**3*(-12*a*b+24*a+12*b+201)
      + j**2*(-3*a**2*b**2+3*a**2*b+6*a**2-3*a*b**2+24*a*b-36*a+6*b**2-27*b-244)
      + j*(-3*a**2*b**2-3*a**2*b+18*a**2-9*a*b**2-15*a*b+66*a+12*b**2+18*b+36))
Q3=K3+pert
print("identity check Q3 == (a-1)(b-2)H - 720C(j,6):", sp.simplify(Q3-Q3id)==0)

def v2(n):
    n=int(n)
    if n==0: return 10**9
    c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n): return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N: return 10**9
    return s2(k)+s2(N-k)-s2(N)

Qf=sp.lambdify((a,b,j),Q3,'math')
def T3(av,bv,jv):
    q0=int(round(Qf(av,bv,0))); qj=int(round(Qf(av,bv,jv)))
    return v2C(av+4,jv)+v2C(bv+3,jv)+v2(qj)-v2(q0)
def B_from_T(av,bv,jv):
    return jv-2*s2(jv)+2*T3(av,bv,jv)   # = Delta(j)

# Empirically: in c2, T(j)>=1-v2(j). Test candidate bounds for c3, split by a parity.
# Look for min of (T3(j)+v2(j)) over each case, and where Delta hits 0.
import collections
for parity,name in [(0,'a EVEN'),(1,'a ODD')]:
    minTplus=collections.defaultdict(lambda:10**9)  # keyed by v2(j)
    tie_js=collections.Counter()
    for mv in range(3,40):
        n=2*mv
        for av in range(1,n+1):
            if av%2!=parity: continue
            for bv in range(3,av+1):
                if n-av-bv!=3: continue
                # find Delta over interior, baseline = min
                deltas={jv:B_from_T(av,bv,jv) for jv in range(1,bv+1)}
                deltas[0]=0
                V=min(deltas.values())
                for jv,d in deltas.items():
                    if d==V and jv>0: tie_js[jv]+=1
                for jv in range(1,bv+1):
                    minTplus[v2(jv)]=min(minTplus[v2(jv)], T3(av,bv,jv)+v2(jv))
    print(f"\n== {name} ==")
    print("  min over cases of T3(j)+v2(j), keyed by v2(j):", dict(minTplus))
    print("  tie j-values (relative to j=0 baseline... note offset):", dict(tie_js))
