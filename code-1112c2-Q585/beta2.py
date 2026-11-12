"""Second tableau<->hw bijection: start from b (shape lam), COLUMN-insert the letters of read(u_mu) from right to
left (column insertion of the reversed word gives the same P as row-inserting the word); record each new box
with the inserted letter.  Candidate bijection hw(B(lam)⊗B(mu))_nu -> LR(nu/lam, mu)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112-Q585')
from commutor import LR, parts, u, rows_read
from flock import is_lr, content
def col_insert(P, x):
    """column insertion: insert x into column 0; bump smallest entry >= x; returns new box (r,c)"""
    cols = {}
    for r,row in enumerate(P):
        for c,v in enumerate(row): cols.setdefault(c,[]).append(v)
    c = 0
    while True:
        col = cols.get(c, [])
        j = next((j for j,y in enumerate(col) if y >= x), None)
        if j is None:
            col.append(x); cols[c] = col; newbox = (len(col)-1, c); break
        col[j], x = x, col[j]; c += 1
    R = max(len(v) for v in cols.values())
    P2 = [[cols[c][r] for c in sorted(cols) if r < len(cols[c])] for r in range(R)]
    return P2, newbox
def beta_col(b, lam, mu):
    P = [list(r) for r in b]; T = {}
    for x in reversed(rows_read(u(mu))):
        P, (r, c) = col_insert(P, x); T[(c, r)] = x
    return T, P
if __name__=='__main__':
    bad=0; tot=0; inj={}
    for a in range(1,6):
      for cc in range(1,6):
        for lam in parts(a,4):
          for mu in parts(cc,4):
            for nu,bs in LR(lam,mu,4).items():
              for b in bs:
                T,P=beta_col(b,lam,mu); tot+=1
                nu_t=tuple(x for x in nu if x)
                okshape=set(T)=={(c,r) for r,l in enumerate(nu_t) for c in range(lam[r] if r<len(lam) else 0,l)}
                if not (is_lr(T) and content(T)==mu and okshape): bad+=1
                inj.setdefault((lam,mu,nu,tuple(sorted(T.items()))),[]).append(b)
    print('beta_col: bad',bad,'of',tot,'; non-injective keys',sum(len(v)>1 for v in inj.values()))
