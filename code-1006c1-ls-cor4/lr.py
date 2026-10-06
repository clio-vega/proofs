"""Littlewood-Richardson coefficients by direct LR-skew-tableau enumeration.
A third, wholly independent mechanism: no Schubert polynomials, no divided
differences, no Bruhat order."""
import itertools
def lr_coeff(lam, mu, nu):
    """c^lam_{mu,nu} = # LR skew tableaux of shape lam/mu and content nu"""
    N=max(len(lam),len(mu))+8
    lam=list(lam)+[0]*(N-len(lam)); mu=list(mu)+[0]*(N-len(mu))
    rows=len([x for x in lam if x>0])
    if sum(lam)!=sum(mu)+sum(nu): return 0
    if any(mu[i]>lam[i] for i in range(len(lam))): return 0
    cells=[(i,j) for i in range(rows) for j in range(mu[i],lam[i])]
    k=len(nu); cnt=0
    def ok_partial(T):
        # semistandard: rows weakly increase, columns strictly increase
        for (i,j),v in T.items():
            if (i,j-1) in T and T[(i,j-1)]>v: return False
            if (i-1,j) in T and T[(i-1,j)]>=v: return False
            if j-1>=mu[i] and (i,j-1) not in T: pass
        return True
    def rec(idx, T, used):
        nonlocal cnt
        if idx==len(cells):
            cnt+=1; return
        (i,j)=cells[idx]
        for v in range(1,k+1):
            if used[v-1]>=nu[v-1]: continue
            if (i,j-1) in T and T[(i,j-1)]>v: continue
            if (i-1,j) in T and T[(i-1,j)]>=v: continue
            T[(i,j)]=v; used[v-1]+=1
            # lattice word condition checked at the end via reverse reading word
            rec(idx+1,T,used)
            del T[(i,j)]; used[v-1]-=1
    # enumerate all SSYT then filter by the lattice (Yamanouchi) condition
    res=[]
    def rec2(idx,T,used):
        if idx==len(cells):
            # reverse reading word: right-to-left, top-to-bottom
            word=[]
            for i in range(rows):
                for j in range(lam[i]-1, mu[i]-1, -1):
                    word.append(T[(i,j)])
            c=[0]*(k+1)
            for v in word:
                c[v]+=1
                if v>1 and c[v]>c[v-1]: return
            res.append(1); return
        (i,j)=cells[idx]
        for v in range(1,k+1):
            if used[v-1]>=nu[v-1]: continue
            if (i,j-1) in T and T[(i,j-1)]>v: continue
            if (i-1,j) in T and T[(i-1,j)]>=v: continue
            T[(i,j)]=v; used[v-1]+=1
            rec2(idx+1,T,used)
            del T[(i,j)]; used[v-1]-=1
    rec2(0,{},[0]*k)
    return len(res)
