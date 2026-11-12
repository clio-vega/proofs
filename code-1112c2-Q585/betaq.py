"""beta_Q: w = read(b)+read(u_mu) (hw, weight nu).  Q(w) standard of shape nu; its entries > |lam| form a
standard skew tableau S on nu/lam, which rectifies to Q(read(u_mu)) (Knuth/jdt).  Relabel entry |lam|+t by
1 + (row of t in Q(read(u_mu))).  Result: LR tableau on nu/lam, content mu (it rectifies to u_mu)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112-Q585')
from commutor import LR, parts, u, rows_read, inv_rsk
from crystal import rsk
from flock import is_lr, content
def beta_q(b, lam, mu):
    wb = rows_read([list(r) for r in b]); m = len(wb)
    U = []; k = 1
    for l in mu: U.append(list(range(k, k+l))); k += l
    wu = inv_rsk(u(mu), U)                      # Knuth-equivalent to u_mu, Q(wu) = row-superstandard U
    P, Q = rsk(wb + wu); _, Qu = rsk(wu)
    rowof = {t: r for r, row in enumerate(Qu) for t in row}
    return {(c, r): rowof[t-m]+1 for r, row in enumerate(Q) for c, t in enumerate(row) if t > m}
if __name__=='__main__':
    bad=0; tot=0; inj={}
    for a in range(1,6):
      for cc in range(1,6):
        for lam in parts(a,4):
          for mu in parts(cc,4):
            for nu,bs in LR(lam,mu,4).items():
              nu_t=tuple(x for x in nu if x)
              for b in bs:
                T=beta_q(b,lam,mu); tot+=1
                okshape=set(T)=={(c,r) for r,l in enumerate(nu_t) for c in range(lam[r] if r<len(lam) else 0,l)}
                if not (is_lr(T) and content(T)==mu and okshape): bad+=1
                inj.setdefault((lam,mu,nu,tuple(sorted(T.items()))),[]).append(b)
    print('beta_Q: bad',bad,'of',tot,'; non-injective keys',sum(len(v)>1 for v in inj.values()))
