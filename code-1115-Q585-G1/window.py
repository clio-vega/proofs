"""Size of the crossing window in B/C diamonds: x's journey in P1 vs P2 (and y's) as (region,type) lists;
window = middle part after removing the longest common prefix and suffix.  Classifies windows by (len_x1,len_x2,len_y1,len_y2)
and by the type-words in the window.  Restricted to the first 40 cells for speed."""
from instr import *
from local import sig, pre, key_of
from collections import Counter
st = Counter(); words = Counter()
def win(a, b):
    p = pre(a, b); s = pre(a[p:][::-1], b[p:][::-1]); return a[p:len(a)-s], b[p:len(b)-s]
ncell = 0
for key, b0 in census():
    ncell += 1
    if ncell > 40: break
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    stack=[b0]; seen=set(); ns=0
    while stack and ns < 60:
        b = stack.pop(); k = key_of(b)
        if k in seen: continue
        seen.add(k); ns += 1; cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
        for x in cs:
            for y in cs:
                if x == y: continue
                b1 = copy.deepcopy(b); _, _, s1x = mig_logged(b1,0,x,DV['B']); _, _, s1y = mig_logged(b1,0,y,DV['C'])
                b2 = copy.deepcopy(b); _, _, s2y = mig_logged(b2,0,y,DV['C']); _, _, s2x = mig_logged(b2,0,x,DV['B'])
                wx = win(sig(s1x), sig(s2x)); wy = win(sig(s1y), sig(s2y))
                if not (wx[0] or wx[1] or wy[0] or wy[1]): continue
                st[(len(wx[0]), len(wx[1]), len(wy[0]), len(wy[1]))] += 1
                words[(''.join(t for _, t in wx[0]), ''.join(t for _, t in wx[1]), ''.join(t for _, t in wy[0]), ''.join(t for _, t in wy[1]))] += 1
print('window lengths (x in P1, x in P2, y in P1, y in P2):'); [print(' ', k, v) for k, v in st.most_common()]
print('type words:'); [print(' ', k, v) for k, v in words.most_common()]
