"""NEGATIVE CONTROLS.  Each must MOVE a number in the table; moved==0 invalidates the control."""
from graph import *
from hookgraph import hooks, conj
from itertools import combinations

def tab(edgefn, lo=1, hi=10):
    out={}
    for n in range(lo,hi+1):
        c=b=0
        for mu in partitions(n):
            V=vertices(mu); E=edgefn(mu)
            comps,comp,colour,odd=components_and_bipartite(V,E)
            if len(comps)==1: c+=1
            if odd: b+=1
        out[n]=(c,b)
    return out

base = tab(lambda mu: edges(mu))
print('BASE          (connected, non-bipartite) per n:', base)

# Control 1: delete all type-(I) edges.  Prediction: connectivity collapses for ell>=2.
c1 = tab(lambda mu: edges(mu, use_I=False))
print('CTRL-1 no (I):', c1, ' moved =', sum(1 for n in base if base[n]!=c1[n]))

# Control 2: delete all type-(II) edges.
c2 = tab(lambda mu: edges(mu, use_II=False))
print('CTRL-2 no(II):', c2, ' moved =', sum(1 for n in base if base[n]!=c2[n]))

# Control 3: replace the EXACT type-(II) condition "b-d1-d2 is a bead" by the paper's
# SUFFICIENT condition "d1+d2 > h_i1".  Prediction: (3,1,1,1,1) loses its triangle.
def edges_suff(mu):
    ell,bs,is_bead = maya(mu)
    E=[]
    holes={b:[u for u in range(-ell,b) if not is_bead(u)] for b in bs}
    for c,b in combinations(sorted(bs),2):
        for bp in holes[c]:
            if not is_bead(b+c-bp):
                E.append(((b,bp),(c,bp),sum(1 for v in range(c+1,b) if is_bead(v)),'I'))
    for b in bs:
        for bp,e in combinations(holes[b],2):
            d1,d2 = b-bp, b-e           # d1 > d2
            if d1+d2 > (b - holes[b][0]):    # h_i1 = b - (smallest hole below b)
                E.append(((b,bp),(b,e),1+sum(1 for v in range(bp+1,e) if is_bead(v)),'II'))
    return E
c3 = tab(edges_suff)
print('CTRL-3 suff :', c3, ' moved =', sum(1 for n in base if base[n]!=c3[n]))
V=vertices((3,1,1,1,1))
print('   (3,1,1,1,1): exact non-bipartite =',bool(components_and_bipartite(V,edges((3,1,1,1,1)))[3]),
      '| sufficient-only non-bipartite =',bool(components_and_bipartite(V,edges_suff((3,1,1,1,1)))[3]))

# Control 4: perturb Lemma 5 -- claim "line of length>=3 always has a triangle".
def L5_wrong(hs): return len(hs)>=3
def L5_right(hs):
    r=len(hs)
    return r>=4 or (r==3 and hs[1]+hs[2]!=hs[0])
def has_tri(hs):
    S=set(hs); return any(hs[q]+hs[s] not in S for q,s in combinations(range(1,len(hs)),2))
w=r=0; wit=[]
for n in range(1,12):
    for mu in partitions(n):
        h=hooks(mu); mp=conj(mu)
        lines=[[h[(i,j)] for j in range(1,mu[i-1]+1)] for i in range(1,len(mu)+1)]
        lines+=[[h[(i,k)] for i in range(1,mp[k-1]+1)] for k in range(1,mu[0]+1)]
        for L in lines:
            if L5_wrong(L)!=has_tri(L):
                w+=1
                if len(wit)<4: wit.append((mu,L))
            if L5_right(L)!=has_tri(L): r+=1
print(f'CTRL-4 Lemma5 without the r=3 exception: {w} wrong predictions (moved={w>0}); correct form: {r}')
print('   witnesses:',wit)
