"""Cor 3.2 bijectivity (injective per boundary) + tie-break lesion (pick the OTHER candidate)."""
import sys
from flock import *
import migrate
neg=lambda v: tuple(-a for a in v)
def unique_flock(b,X,E,N):
    pos=b.nest[X]['pos']; bx={box_of(b,X,b.tiles[t][1],E,N):t for t in pos.values()}
    sh=tuple(sorted(Counter(r for _,r in bx).values(),reverse=True))
    L=lr_tableaux(set(bx),sh); assert len(L)==1
    return {bx[k]:v for k,v in L[0].items()}
def run(n):
    seen={}; ok=fail=0
    for lab in puzzles(n,{}):
        left,right,bottom=side_edges(n); k1=sum(lab[E(*e)] for e in left)
        if k1 in (0,n): continue
        bnd=tuple(tuple(lab[E(*e)] for e in s) for s in (left,right,bottom))
        b=Board(build_mosaic(n,lab,SIGN,RULE)); A=b.nest[0]
        E_,N_=neg(A['u']),neg(A['w'])
        try:
            fl=unique_flock(b,0,E_,N_) if A['pos'] else {}
            boxes={} if not fl else migrate_flock(b,0,fl,E_,N_,'S')[3]
            ok+=1
        except Exception as ex: fail+=1; continue
        key=(bnd,tuple(sorted(boxes.items())))
        seen[key]=seen.get(key,0)+1
    coll=sum(v-1 for v in seen.values())
    return ok,fail,coll
n=int(sys.argv[1])
print('n',n,'(ok, fail, collisions within a boundary)',run(n))
if len(sys.argv)>2:
    # lesion: invert the tie-break
    src=open('migrate.py').read()
    src=src.replace("cands = [zz for zz in cands if across & zz[2]]","cands = [zz for zz in cands if not (across & zz[2])]")
    exec(compile(src,'migrate_lesion','exec'),migrate.__dict__)
    import flock; flock.Board=migrate.Board
    globals()['Board']=migrate.Board
    print('LESION (other tie choice): (ok, fail, collisions)',run(n))
