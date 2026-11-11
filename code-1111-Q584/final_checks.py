import itertools, sys
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import D
from star import UP, DN, X
# (a) 4 unit lattice vectors summing to 0: always two opposite pairs; simple cyclic orders = rhombi
V=[D[i] for i in range(6)]
bad=0; tot=0; simple=set()
for w in itertools.product(range(6),repeat=4):
    if sum(D[i][0] for i in w)==0 and sum(D[i][1] for i in w)==0:
        tot+=1
        ms=sorted(w)
        pairs = all(((i+3)%6) in w for i in w) and sorted(w)==sorted([(i+3)%6 for i in w])
        if not pairs: bad+=1
        backtrack=any((w[k]-w[(k+1)%4])%6==3 for k in range(4))
        if not backtrack: simple.add(tuple(sorted({w[0]%3,w[1]%3})))
print('(a) zero-sum 4-step words:',tot,' not two opposite pairs:',bad,' non-backtracking class pairs:',sorted(simple))
# (b) sigma_W o C on unit triangles, sigma_W: dir d -> (2-d)%6 ; C: 0<->1, X fixed.
# up triangle (dir0,dir2,dir4) labels (l0,l2,l4) -> image triangle has edges in dirs 2,0,4 reversed orientation -> a DOWN-class? compute directly
def img_tri(kind,labs):
    dirs=(0,2,4) if kind=='U' else (1,3,5)
    e=[((2-d)%6, {0:1,1:0,X:X}[l]) for d,l in zip(dirs,labs)]
    # reflected triangle traversed cw -> reverse & flip dir by 3 to get ccw
    e=[((d+3)%6,l) for d,l in reversed(e)]
    # rotate list to start at the smallest dir and identify orientation
    e=sorted(e); ds=tuple(d for d,_ in e)
    return ('U' if ds==(0,2,4) else 'D' if ds==(1,3,5) else '?', tuple(l for _,l in e))
ok=True
for kind,S,T in (('U',UP,DN),('D',DN,UP)):
    for t in list(S)+[(X,X,X)]:
        k2,l2=img_tri(kind,t)
        inside = l2 in (T|{(X,X,X)})
        ok &= inside
        if t==(X,X,X): print('(b) K-piece',kind,'->',k2,l2)
print('(b) sigma_W o C maps {UP+K} <-> {DN+K} bijectively:',ok)
