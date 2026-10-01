from regions import *

print("=== INSTRUMENT CHECK 1: reproduce the KNOWN FAILURE of (*) between slice members ===")
print("    n=6, mu=(0,3), lam=(4,5), slice S_nu=5, nu=(0,5) vs rho=(2,3)")
n, mu, lam = 6, (0,3), (4,5)
u = 5; b = u - sum(mu)
nn, A, B, uu, Tm, Tp = cyl_params(n, mu, lam, b)
print(f"    A={A} B={B} C={B+1} D={A+n-1} u={uu} T=[{Tm},{Tp}] tau1={A} tau2={uu-1-B}")
members = {}
for t in range(Tm, Tp+1):
    at,bt,ct,dt = endpoints(t,n,A,B,uu)
    members[t] = conv_interval(at,bt,ct,dt)
    print(f"    t={t}: [a,b]=[{at},{bt}] [c,d]=[{ct},{dt}] region={region_of(t,A,uu-1-B)}  trap={as_seq(members[t])}")
bad = star_fails(members[0], members[2])
print(f"    (*) for t=0 vs t=2 : failures = {bad}")
assert bad, "INSTRUMENT BROKEN: cannot reproduce the known failure"
print("    OK: known failure reproduced.")

print()
print("=== INSTRUMENT CHECK 2: G_R individually PF2 on this slice ===")
G, tset = region_sums(n, A, B, uu, Tm, Tp)
for R in ['I','II','III','IV']:
    print(f"    {R}: t={tset[R]} G={as_seq(G[R])} PF2={is_pf2(G[R])}")
assert all(is_pf2(G[R]) for R in G), "INSTRUMENT BROKEN: a G_R is not PF2"
print("    OK: all four PF2.")

print()
print("=== INSTRUMENT CHECK 3: G = sum G_R equals brute-force triple count ===")
# brute force: count (t, kappa1, kappa2) with constraints, graded by kappa1+kappa2
tot = {}
for t in range(Tm, Tp+1):
    for s,v in conv_interval(*endpoints(t,n,A,B,uu)).items():
        tot[s] = tot.get(s,0)+v
summed = {}
for R in G:
    for s,v in G[R].items(): summed[s]=summed.get(s,0)+v
print(f"    total={as_seq(tot)}  sum of regions={as_seq(summed)}  equal={tot==summed}  PF2(G)={is_pf2(tot)}")
assert tot == summed
print("    OK.")

print()
print("=== INSTRUMENT CHECK 4: regions II/IV constant-support claim ===")
# region II: a_t+c_t=u, b_t+d_t=B+D ; region IV: a_t+c_t=A+C, b_t+d_t=u+n-2
fails=0; checked=0
for n_ in range(2,9):
  for A_ in range(-6,7):
    for B_ in range(A_,A_+9):
      C_,D_ = B_+1, A_+n_-1
      for u_ in range(-4,20):
        for t in range(-8,20):
          at,bt,ct,dt = endpoints(t,n_,A_,B_,u_)
          R = region_of(t,A_,u_-1-B_)
          if at>bt or ct>dt: continue
          checked+=1
          if R=='II' and not (at+ct==u_ and bt+dt==B_+D_): fails+=1
          if R=='IV' and not (at+ct==A_+C_ and bt+dt==u_+n_-2): fails+=1
print(f"    checked={checked} fails={fails}")
