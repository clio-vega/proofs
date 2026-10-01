from regions import *
from itertools import combinations
ORD=['I','II','III','IV']
def tot_of(G):
    tot={}
    for R in G:
        for s,v in G[R].items(): tot[s]=tot.get(s,0)+v
    return tot

print("=== classify all (G) failures, cylindric n<=11, d<=14 ===")
fails=[]; nslice=0; npair=0; A_fail=0; A_fail_wit=[]
nregion_hist={}
twoRegion_pairfail=0
for n in range(2,12):
  for mu1 in range(0,n+1):
    for mu2 in range(mu1, mu1+n+1):
      for lam1 in range(mu1, mu1+16):
        for lam2 in range(max(mu2,lam1), lam1+n+1):
          d=(lam1-mu1)+(lam2-mu2)
          if d>14 or d<1: continue
          if not valid_shape(n,(mu1,mu2),(lam1,lam2)): continue
          for b in range(0,d+1):
            nn,A,B,u,Tm,Tp=cyl_params(n,(mu1,mu2),(lam1,lam2),b)
            if Tm>Tp: continue
            nslice+=1
            G,tset=region_sums(n,A,B,u,Tm,Tp)
            nz=[R for R in ORD if G[R]]
            nregion_hist[len(nz)]=nregion_hist.get(len(nz),0)+1
            tot=tot_of(G)
            if not is_pf2(tot):
                A_fail+=1
                if len(A_fail_wit)<3: A_fail_wit.append((n,(mu1,mu2),(lam1,lam2),b,as_seq(tot)))
            for R,Rp in combinations(nz,2):
                npair+=1
                bad=star_fails(G[R],G[Rp])
                if bad:
                    fails.append(dict(n=n,mu=(mu1,mu2),lam=(lam1,lam2),b=b,pair=(R,Rp),
                                      bad=bad,gR=as_seq(G[R]),gRp=as_seq(G[Rp]),nz=tuple(nz),
                                      totPF2=is_pf2(tot)))
                    if len(nz)==2: twoRegion_pairfail+=1
print(f"slices={nslice} pairs={npair} (G)-failures={len(fails)}")
print(f"regions-nonempty histogram: {nregion_hist}")
print(f"CONDITION (A) failures (G not PF2): {A_fail}  {A_fail_wit}")
print(f"(G)-failures occurring when only TWO regions nonempty: {twoRegion_pairfail}")
print()
from collections import Counter
print("failures by pair:", Counter(f['pair'] for f in fails))
print("failures by nonempty-region-set:", Counter(f['nz'] for f in fails))
print("in every failure, is the total G still PF2? ", all(f['totPF2'] for f in fails))
print()
print("--- is the 'flat partner' mechanism universal? ---")
# mechanism: the partner sequence is constant on an interval containing supp of the other
mech=0
for f in fails:
    (R,Rp)=f['pair']
    g1=f['gR']; g2=f['gRp']
    # test both orientations: one is constant >0 on a window covering the other's support
    def flat_cover(ga,gb):
        loa,sa=ga; lob,sb=gb
        # gb constant on [loa, loa+len(sa)-1]?
        vals=[sb[k-lob] if 0<=k-lob<len(sb) else 0 for k in range(loa,loa+len(sa))]
        return len(set(vals))==1 and vals[0]>0
    if flat_cover(g1,g2) or flat_cover(g2,g1): mech+=1
print(f"  failures explained by 'one region flat across the other''s support': {mech}/{len(fails)}")
print()
print("--- minimal witnesses (smallest n, then d, then |mu|) ---")
fails.sort(key=lambda f:(f['n'],(f['lam'][0]-f['mu'][0])+(f['lam'][1]-f['mu'][1]),sum(f['mu'])))
for f in fails[:8]:
    d=(f['lam'][0]-f['mu'][0])+(f['lam'][1]-f['mu'][1])
    print(f"  n={f['n']} mu={f['mu']} lam={f['lam']} d={d} b={f['b']} pair={f['pair']} nz={f['nz']}")
    print(f"     g_{f['pair'][0]}={f['gR']}  g_{f['pair'][1]}={f['gRp']}  bad s=(s,lhs,rhs)={f['bad']}  G PF2={f['totPF2']}")
