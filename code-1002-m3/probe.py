"""Structural probes for the aggregate layer profile beta.
P-a  is sum_r beta(r) x^r real-rooted?   (would open an interlacing proof)
P-b  is the set of half-widths Lambda(nu) on a slice gap-free? (generalises the tent lemma)
P-c  pairwise (*) among the LAYER profiles beta_nu: census of failures
P-d  is beta unimodal?
"""
import gen, layers, sys
from collections import Counter
import numpy as np

def realrooted(co):
    co=[float(c) for c in co]
    while co and co[-1]==0: co.pop()
    while co and co[0]==0: co.pop(0)
    if len(co)<=1: return True
    r=np.roots(list(reversed(co)))
    return bool(np.all(np.abs(r.imag) < 1e-7*np.maximum(1.0,np.abs(r.real))))

def star(g,h):
    N=max(len(g),len(h))+2
    G=lambda i: g[i] if 0<=i<len(g) else 0
    H=lambda i: h[i] if 0<=i<len(h) else 0
    for s in range(-1,N):
        if 2*G(s)*H(s) < G(s-1)*H(s+1)+G(s+1)*H(s-1): return False
    return True

st=Counter(); ex=Counter()
wit_rr=[]; wit_gap=[]
for m in (2,3,4,5):
    for n in range(m,10):
        for (mu,lam) in gen.pairs(n,m,11):
            for b in range(0,sum(lam)-sum(mu)+1):
                nus=gen.slice_nus(mu,lam,n,m,b)
                fs=[gen.f_nu(nu,lam,n,m) for nu in nus]
                fs=[f for f in fs if f[1]]
                if not fs: continue
                off,co=gen.slice_sum(mu,lam,n,m,b)
                bet=layers.beta_of(layers.radial(co))
                st[f'm{m}_slices']+=1
                # P-a
                if realrooted(bet): st[f'm{m}_beta_realrooted']+=1
                elif len(wit_rr)<4: wit_rr.append((m,n,mu,lam,b,co,bet))
                if realrooted(co): st[f'm{m}_G_realrooted']+=1
                # P-b  half-widths
                lams=sorted({(len(c)-1)/2 for (o,c) in fs})
                if all(lams[i+1]-lams[i] in (0.5,1.0) for i in range(len(lams)-1)):
                    st[f'm{m}_Lambda_gapfree']+=1
                else:
                    if len(wit_gap)<4: wit_gap.append((m,n,mu,lam,b,lams))
                # P-c pairwise (*) on layer profiles
                bs=[layers.beta_of(layers.radial(c)) for (o,c) in fs]
                okpairs=badpairs=0
                for i in range(len(bs)):
                    for j in range(i+1,len(bs)):
                        if star(bs[i],bs[j]): okpairs+=1
                        else: badpairs+=1
                st[f'm{m}_layerpair_ok']+=okpairs; st[f'm{m}_layerpair_BAD']+=badpairs
                # P-d unimodal
                inc=[bet[i+1]-bet[i] for i in range(len(bet)-1)]
                signs=[1 if x>0 else (-1 if x<0 else 0) for x in inc]
                sg=[s for s in signs if s]
                st[f'm{m}_beta_unimodal' if all(sg[i]>=sg[i+1] for i in range(len(sg)-1)) else f'm{m}_beta_NOTunimodal']+=1
for k in sorted(st): print(f"{k:28s} {st[k]}")
print("\nbeta not real-rooted e.g.:", wit_rr[:2])
print("\nLambda has a gap e.g.:", wit_gap[:2])
