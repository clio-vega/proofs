import sys; sys.path.insert(0,'.')
from commutor import *
def C2_lesion(b,lam,mu,n):   # lesion: drop the outer xi (keep P)
    xb=rows_read(evac([list(r) for r in b],n)); xc=rows_read(evac(u(mu),n))
    x=xc+xb
    return split_rows(x[:sum(mu)],mu), split_rows(x[sum(mu):],lam)
def F_lex(b,lam,mu,nu,n,L1,L2):  # control 1: lex-rank matching
    return sorted(L2[nu])[sorted(L1[nu]).index(b)]
def F_rev(b,lam,mu,nu,n,L1,L2):  # control 2: reverse-lex matching
    return sorted(L2[nu],reverse=True)[sorted(L1[nu]).index(b)]
n=4; les=0; tot=0; cells2=0; dlex=0; drev=0; fixed=0; sym=0
for a in range(1,6):
  for bs_ in range(1,6):
    for lam in parts(a,n):
      for mu in parts(bs_,n):
        L1=LR(lam,mu,n); L2=LR(mu,lam,n)
        for nu,bs in L1.items():
          for b in bs:
            tot+=1; les+= C2_lesion(b,lam,mu,n)!=C1(b,lam,mu,n)
          if len(bs)<2: continue
          cells2+=1
          s={b:C1(b,lam,mu,n)[0] for b in bs}
          dlex+= any(s[b]!=F_lex(b,lam,mu,nu,n,L1,L2) for b in bs)
          drev+= any(s[b]!=F_rev(b,lam,mu,nu,n,L1,L2) for b in bs)
          if lam==mu: sym+=1; fixed+= all(s[b]==b for b in bs)
print('lesion (no outer xi) differs from C1 on',les,'of',tot)
print('c>=2 cells',cells2,'; C differs from lex-matching on',dlex,'; from revlex on',drev)
print('lam==mu cells with c>=2:',sym,'; of these C is identity on',fixed)
