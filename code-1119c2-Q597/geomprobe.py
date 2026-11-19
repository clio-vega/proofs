import sys
sys.argv=[sys.argv[0],'5','1','1']
exec(open('/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py').read().split('\ndef work(w):')[0])
from collections import Counter
st=Counter(); n=0
for key,b in census_sized():
    for bb,DV in states(b,cap=CAP):
        M0=bare(bb)
        for t in M0:
            cl=classes(t)
            cp=tuple(sorted('L0' if c[2:]==(0,0) else 'L1' for c in cl))
            st[(t[0],len(t[1]),cp)]+=1
        # lozenge pairing: each tri shares an edge with how many tris?
        E=index(M0)
        for t in M0:
            if t[0]!='tri': continue
            k=sum(1 for e in tedges(t) for s in E[e] if s!=t and s[0]=='tri')
            st[('tri-tri neighbours',k)]+=1
        n+=1
        if n>=40: break
    if n>=40: break
for k,v in sorted(st.items(),key=str): print(v,k)
