"""Q584 step 0: D6 x {1,C} acting on Q569 tiles. C: 0<->1, 2 fixed."""
import sys; sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import TILES, D
from symm import act, canon, GROUP, NAME, INV, KTW
def C(steps): return [(d,{0:1,1:0,2:2}[l]) for d,l in steps]
def img(nm,g,c):
    s=act(TILES[nm],g); s=C(s) if c else s
    return INV.get(canon(s))
# positive control on orientation from BASIS IMAGES: det of linear map on e0,e1
def lin(g):
    t,p=g
    f=(lambda d:(d+p)%6) if t=='r' else (lambda d:(p-d)%6)  # direction map on lines
    e0,e1=D[f(0)],D[f(1)]
    return e0[0]*e1[1]-e0[1]*e1[0]
ktw=set(KTW)
for c in (0,1):
  for g in GROUP:
    pres = all(img(t,g,c) in ktw for t in KTW)
    perm = {t:img(t,g,c) for t in KTW}
    print(('C∘' if c else '  ')+NAME[g].ljust(8),'det',lin(g),'KTW-preserved',pres,
          '| hex->',img('hex',g,c),' trap->',img('trap',g,c),' K->',img('K',g,c),' tri->',img('tri',g,c),
          ' hex2->',img('hex2',g,c))
