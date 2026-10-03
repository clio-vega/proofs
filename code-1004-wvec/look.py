"""Stay with the object: print real width-vector families, grouped by half-width, and
look at what structure they have.  Also: the y-set Y (box cap hyperplane, M-convex by
A-M-profile-box-slice) and the map y -> w on it."""
import gen, tent
from collections import defaultdict

def data(mu,lam,n,m,b):
    target=sum(mu)+b; rows=[]
    for nu in gen.box(mu,n,m):
        if sum(nu)!=target: continue
        lr=gen.LR(nu,lam,n,m); w=tuple(R-L+1 for (L,R) in lr)
        if any(x<=0 for x in w): continue
        y=tuple(nu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m))
        rows.append((y,w))
    return rows

shown=0
for m in (3,4):
  for n in range(m+3,m+7):
    for (mu,lam) in gen.pairs(n,m,3*(n-m)):
        d=sum(lam)-sum(mu)
        for b in range(0,d+1):
            rows=data(mu,lam,n,m,b)
            if len(rows)<5: continue
            hw={ (sum(w)-m)//2 for (y,w) in rows}
            if len(hw)<3: continue
            off,co=gen.slice_sum(mu,lam,n,m,b)
            g=tuple(tent.lam_at(lam,i,n,m)-tent.lam_at(lam,i-1,n,m)-1 for i in range(m))
            print(f"--- m={m} n={n} mu={mu} lam={lam} b={b}  G={n-m} g={g} sigma={sum(rows[0][0])} "
                  f"|W|={len({w for _,w in rows})} nsummands={len(rows)} PF2={gen.is_pf2(co)}")
            byhw=defaultdict(list)
            for y,w in rows: byhw[(sum(w)-m)//2].append((y,w))
            for L in sorted(byhw):
                print(f"    Lambda={L}: "+"  ".join(f"y{list(y)}->w{list(w)}" for y,w in sorted(byhw[L])))
            print(f"    slice sum = {co}")
            shown+=1
            break
        if shown>=6: break
    if shown>=6: break
  if shown>=6: break
