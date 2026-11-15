"""Draw the first interacting B/C diamond in the census (smallest cell first): start state, x->B;y->C, y->C;x->B,
with the hexagon steps of each journey listed.  Output: diamond-{0,1,2}.png + printed step words."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585')
from instr import *
from local import sig
from draw import draw
for key, b in census():
    A = b.nest[0]; DV = {'B': A['w'], 'C': A['u']}; cs = corners(set(A['pos'])); done = False
    for x in cs:
        for y in cs:
            if x == y: continue
            b1 = copy.deepcopy(b); _, lx1, s1x = mig_logged(b1,0,x,DV['B']); _, ly1, s1y = mig_logged(b1,0,y,DV['C'])
            b2 = copy.deepcopy(b); _, ly2, s2y = mig_logged(b2,0,y,DV['C']); _, lx2, s2x = mig_logged(b2,0,x,DV['B'])
            if sig(s1x) != sig(s2x):
                print(key, 'x', x, 'y', y, 'landings', lx1, ly1, '|', lx2, ly2)
                print(' x->B first :', ''.join(s['typ'] for s in s1x), ' then y->C:', ''.join(s['typ'] for s in s1y))
                print(' y->C first :', ''.join(s['typ'] for s in s2y), ' then x->B:', ''.join(s['typ'] for s in s2x))
                draw(b, 'diamond-0.png'); draw(b1, 'diamond-1.png'); draw(b2, 'diamond-2.png'); done = True; break
        if done: break
    if done: break
