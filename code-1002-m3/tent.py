"""The general-m tent lemma and its two ingredients, checked against the enumerator."""
import gen, sys
from collections import Counter

def lam_at(lam, j, n, m):
    """lam_j for any integer j, using lam_{j+m} = lam_j + n."""
    q, r = divmod(j, m)
    return lam[r] + q*n

def Lambda_formula(nu, lam, n, m):
    sig = (sum(nu) + sum(lam)) / 2
    return sig - sum(max(nu[i], lam_at(lam, i-1, n, m) + 1) for i in range(m))

def eff_box(mu, lam, n, m):
    """predicted: L_i<=R_i for all i  <=>  nu_i <= lam_i  and  nu_i >= lam_{i-2}+2."""
    return [(max(mu[i], lam_at(lam, i-2, n, m) + 2),
             min(gen.nu_next(mu, i, n, m) - 1, lam[i])) for i in range(m)]

st=Counter(); wit=[]
for m in range(2,6):
    for n in range(m, m+5):
        for (mu,lam) in gen.pairs(n,m,8):
            rng = eff_box(mu,lam,n,m)
            for nu in gen.box(mu,n,m):
                off,co = gen.f_nu(nu,lam,n,m); nz = bool(co); st['nu']+=1
                if nz:
                    if abs(Lambda_formula(nu,lam,n,m) - (len(co)-1)/2) < 1e-9: st['T1_ok']+=1
                    else:
                        st['T1_BAD']+=1
                        if len(wit)<3: wit.append(('T1',m,n,mu,lam,nu,Lambda_formula(nu,lam,n,m),(len(co)-1)/2))
                inbox = all(rng[i][0] <= nu[i] <= rng[i][1] for i in range(m))
                if inbox==nz: st['T2_ok']+=1
                else:
                    st['T2_BAD']+=1
                    if len(wit)<6: wit.append(('T2',m,n,mu,lam,nu,inbox,nz,rng))
print(dict(st)); [print('  ',w) for w in wit[:4]]

print("\n-- refusal controls (perturbed formulas must be REFUSED) --")
for name, pert in [('T1 with lam_{i-1}+2', 2), ('T1 with lam_{i-1}+0', 0)]:
    bad=tot=0
    for m in (3,4):
        for n in range(m,m+3):
            for (mu,lam) in gen.pairs(n,m,7):
                for nu in gen.box(mu,n,m):
                    off,co=gen.f_nu(nu,lam,n,m)
                    if not co: continue
                    tot+=1
                    v=(sum(nu)+sum(lam))/2 - sum(max(nu[i], lam_at(lam,i-1,n,m)+pert) for i in range(m))
                    if abs(v-(len(co)-1)/2)>1e-9: bad+=1
    print(f"   {name}: refused {bad}/{tot}")
for name, pert in [('T2 with lam_{i-2}+3',3), ('T2 with lam_{i-2}+1',1)]:
    bad=tot=0
    for m in (3,4):
        for n in range(m,m+3):
            for (mu,lam) in gen.pairs(n,m,7):
                rng=[(max(mu[i], lam_at(lam,i-2,n,m)+pert), min(gen.nu_next(mu,i,n,m)-1, lam[i])) for i in range(m)]
                for nu in gen.box(mu,n,m):
                    tot+=1
                    nz=bool(gen.f_nu(nu,lam,n,m)[1])
                    if all(rng[i][0]<=nu[i]<=rng[i][1] for i in range(m)) != nz: bad+=1
    print(f"   {name}: refused {bad}/{tot}")
