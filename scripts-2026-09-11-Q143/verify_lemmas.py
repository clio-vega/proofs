from graph import *
from hookgraph import hooks, conj, hook_edges
from itertools import combinations

# ---- Lemma 2 (hubs) ----
fail=0
for n in range(1,11):
    for mu in partitions(n):
        E={frozenset((x,y)) for (x,y,ty) in hook_edges(mu)}; mp=conj(mu)
        for i in range(1,len(mu)+1):
            for m in range(2,mu[i-1]+1):
                if frozenset(((i,1),(i,m))) not in E: fail+=1; print('hub row fail',mu,i,m)
        for k in range(1,(mu[0] if mu else 0)+1):
            for ip in range(2,mp[k-1]+1):
                if frozenset(((1,k),(ip,k))) not in E: fail+=1; print('hub col fail',mu,k,ip)
print('Lemma 2 (hubs): failures =',fail,'  over all mu, n<=10')

# ---- Lemma 5 (line triangle) : abstract statement on the hook sequences ----
def line_has_triangle(hs):
    """hs strictly decreasing. triangle {1,q,s} exists iff some 2<=q<s with hs[q]+hs[s] not in hs"""
    S=set(hs); r=len(hs)
    return any(hs[q]+hs[s] not in S for q,s in combinations(range(1,r),2))
def lemma5_predicts(hs):
    r=len(hs)
    if r<3: return False
    if r==3: return hs[1]+hs[2]!=hs[0]
    return True
f2=0; tot=0
for n in range(1,12):
    for mu in partitions(n):
        h=hooks(mu); mp=conj(mu)
        for i in range(1,len(mu)+1):
            hs=[h[(i,j)] for j in range(1,mu[i-1]+1)]; tot+=1
            if line_has_triangle(hs)!=lemma5_predicts(hs): f2+=1; print('L5 row fail',mu,i,hs)
        for k in range(1,(mu[0] if mu else 0)+1):
            hs=[h[(i,k)] for i in range(1,mp[k-1]+1)]; tot+=1
            if line_has_triangle(hs)!=lemma5_predicts(hs): f2+=1; print('L5 col fail',mu,k,hs)
print(f'Lemma 5 (line triangle): failures = {f2} over {tot} lines (rows+columns), n<=11')

# ---- Theorem 6: n>=6 => triangle; which line supplies it ----
print()
print(' n  #mu  has-triangle  row1-or-col1-suffices   residual mu needing another line')
for n in range(1,12):
    ps=list(partitions(n)); tri=0; r1c1=0; resid=[]
    for mu in ps:
        h=hooks(mu); mp=conj(mu)
        lines=[[h[(i,j)] for j in range(1,mu[i-1]+1)] for i in range(1,len(mu)+1)]
        lines+=[[h[(i,k)] for i in range(1,mp[k-1]+1)] for k in range(1,mu[0]+1)]
        any_t=any(line_has_triangle(L) for L in lines)
        row1=line_has_triangle([h[(1,j)] for j in range(1,mu[0]+1)])
        col1=line_has_triangle([h[(i,1)] for i in range(1,len(mu)+1)])
        if any_t: tri+=1
        if row1 or col1: r1c1+=1
        elif any_t: resid.append(mu)
        elif n>=6: print('   *** NO TRIANGLE AT ALL:',mu)
    print(f'{n:2d} {len(ps):4d} {tri:9d} {r1c1:14d}        {resid}')
