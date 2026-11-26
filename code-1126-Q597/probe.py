import sys
sys.argv=['x','5','1','24']
src=open('/home/clio/projects/proofs/code-1125c2-Q597/bsall.py').read().split("\ndef work(w):")[0]
exec(src)
from collections import Counter
n=0
for idx,(key,b) in enumerate(census_sized()):
    if idx!=3: continue
    for b,DV in states(b,cap=CAP):
        cs=corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x==y: continue
                b1,Jx=mvl(b,x,'B',DV); b2,Jy=mvl(b,y,'C',DV)
                b12,Jyx=mvl(b1,y,'C',DV); b21,Jxy=mvl(b2,x,'B',DV)
                if Jx==Jxy and Jy==Jyx: continue
                print('x',x,'y',y)
                for nm,A,Bq in (('x',Jx,Jxy),('y',Jy,Jyx)):
                    mA,mB,p,q=middles(A,Bq)
                    print(nm,'p',p,'len',len(A),len(Bq),'A',ws(typed(A)),'| Ar',ws(typed(Bq)),'mid',len(mA),len(mB))
                print(repr(Jx[0])[:300])
                n+=1
                if n>4: sys.exit()
