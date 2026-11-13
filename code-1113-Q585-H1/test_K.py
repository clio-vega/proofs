"""Hypothesis K1: kappa(P) = hw( revrot(Tao(P)) ), revrot = rotate 180 deg and v -> l+1-v, hw = raise by e_i to highest weight
(reading word: rows top->bottom, each left->right; e_i flips leftmost unmatched i+1, matching (i+1 ... i))."""
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts

def word_cells(T):
    rows=sorted({r for c,r in T},reverse=True)
    return [(c,r) for r in rows for c in sorted(c for c,rr in T if rr==r)]
def e_i(T,i):
    cells=word_cells(T); stack=[]; unmatched=[]
    for x in cells:
        v=T[x]
        if v==i+1: stack.append(x)
        elif v==i and stack: stack.pop()
    if not stack: return None
    T=dict(T); T[stack[0]]=i; return T
def hw(T):
    m=max(T.values()) if T else 0; changed=True
    while changed:
        changed=False
        for i in range(1,m):
            while True:
                T2=e_i(T,i)
                if T2 is None: break
                T=T2; changed=True
    return T
def revrot(T,l):
    return {(-c,-r):l+1-v for (c,r),v in T.items()}
st={}
def bump(k): st[k]=st.get(k,0)+1
for a in range(1,5):
 for c_ in range(1,5):
  for lam in parts(a,4):
   for mu in parts(c_,4):
    for nu,bs in LR(lam,mu,4).items():
      if len(bs)<2: continue
      nt=trim(nu)
      for bb in bs:
        T=beta_cm(bb,lam,mu,nt)
        k0=max(nt[0],1)+1; k1=len(nt)+1
        b=canonical(nt,k0,k1); B=b.nest[1]; EB,NB=B['u'],B['w']
        F={B['pos'][(c,r)]:v for (c,r),v in T.items()}
        migrate_flock(b,1,F,EB,NB,'E')
        _,tao=readout(b,0,1,'N'); bump('n')
        l=len(lam)
        if norm(hw(revrot(tao,l)))==norm(T): bump('K1')
        if norm(revrot(tao,l))==norm(T): bump('revrot alone')
        # control: lex-like wrong choice: other LR tableau of same shape
        if is_lr(hw(revrot(tao,l))): bump('K1 output is LR')
print(st)
