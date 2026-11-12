from flock import *
import copy, sys
from collections import Counter
neg=lambda v: tuple(-a for a in v)
def unique_flock(b,X,E,N):
    pos=b.nest[X]['pos']
    bx={box_of(b,X,b.tiles[t][1],E,N):t for t in pos.values()}
    sh=tuple(sorted(Counter(r for _,r in bx).values(),reverse=True))
    L=lr_tableaux(set(bx),sh); assert len(L)==1
    return {bx[k]:v for k,v in L[0].items()} and {bx[k]:v for k,v in L[0].items()}
n=int(sys.argv[1]); res=Counter()
for lab in puzzles(n,{}):
    left,right,bottom=side_edges(n)
    k1=sum(lab[E(*e)] for e in left)
    if k1 in (0,n): continue
    for oname in ('E,N','-E,-N','N,E','-N,-E'):
      for dname in 'NESW':
        b=Board(build_mosaic(n,lab,SIGN,RULE)); A=b.nest[0]
        if not A['pos']: continue
        u,w=A['u'],A['w']
        E_,N_={'E,N':(u,w),'-E,-N':(neg(u),neg(w)),'N,E':(w,u),'-N,-E':(neg(w),neg(u))}[oname]
        # consistency: opens NE in frame iff (E,N) in {(u,w),(w,u)} -> direction N/E
        if (oname in ('E,N','N,E')) != (dname in 'NE'): continue
        try:
            fl=unique_flock(b,0,E_,N_)
            dst,out,fr2,boxes,st=migrate_flock(b,0,fl,E_,N_,dname)
            res[(oname,dname,'ok->'+'ABC'[dst],is_lr(boxes))]+=1
        except Exception as ex:
            res[(oname,dname,'FAIL',str(ex)[:40])]+=1
for k,v in sorted(res.items()): print(k,v)
