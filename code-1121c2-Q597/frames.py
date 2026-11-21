"""PROVE 11-21 c2: the nest frames (u, w) of nests A, B, C over the <=5 census (FIN: a packed nest-B/C rhombus never has
the nest-A class, i.e. its Lambda_0 edge is not +-(0,1)).   usage: python3 frames.py 5 1 1"""
import os
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
st = Counter()
for idx, (key, b) in enumerate(census_sized()):
    for i, info in enumerate(b.nest):
        u, w = tuple(str(x) for x in info['u']), tuple(str(x) for x in info['w'])
        st[('nest', 'ABC'[i], 'u', u, 'w', w)] += 1
for k in sorted(st, key=str): print(st[k], k)
