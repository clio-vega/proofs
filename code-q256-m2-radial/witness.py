from regions import *
n,mu,lam,b = 10,(0,5),(5,12),4
print("shape valid:",valid_shape(n,mu,lam), " d=",(lam[0]-mu[0])+(lam[1]-mu[1]))
nn,A,B,u,Tm,Tp = cyl_params(n,mu,lam,b)
C,D = B+1, A+n-1
tau1, tau2 = A, u-1-B
print(f"n={n} mu={mu} lam={lam} b={b}  A={A} B={B} C={C} D={D} u={u} T=[{Tm},{Tp}] tau1={tau1} tau2={tau2}")
G,tset = region_sums(n,A,B,u,Tm,Tp)
for t in range(Tm,Tp+1):
    at,bt,ct,dt = endpoints(t,n,A,B,u)
    print(f"  t={t}: [{at},{bt}]x[{ct},{dt}] region={region_of(t,tau1,tau2)} trap={as_seq(conv_interval(at,bt,ct,dt))}")
print()
for R in ['I','II','III','IV']:
    print(f"  G_{R}: t={tset[R]}  seq={as_seq(G[R])}  PF2={is_pf2(G[R])}")
tot={}
for R in G:
    for s,v in G[R].items(): tot[s]=tot.get(s,0)+v
print()
print("  G total =", as_seq(tot), " PF2(G) =", is_pf2(tot))
print("  concave(G_I)?", end=" ")
g=G['I']; lo,hi=min(g),max(g)
print([(s,2*get(g,s)-get(g,s-1)-get(g,s+1)) for s in range(lo,hi+1)])
print()
print("  (*) for I-III failures:", star_fails(G['I'],G['III']))
print()
print("  ** KEY: is G_III constant on the support of G_I? **")
print("     supp G_I =",(min(G['I']),max(G['I'])), " supp G_III =",(min(G['III']),max(G['III'])))
print("     G_III values:",[get(G['III'],s) for s in range(min(G['I']),max(G['I'])+1)])
