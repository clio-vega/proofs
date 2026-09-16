"""For each admissible support Z, the value solutions form a torsor: sigma = 1_Z * chi
with chi a character of Z^Z / Lattice(relations).  Compute the lattice's Smith form."""
import sys; sys.path.insert(0,'.')
from twobead import words, relations_2Bstar
from supports import enumerate_supports
from sympy import Matrix
from analyse import Cset

def lattice_invariants(Z, d, ws):
    Zl = sorted(Z); idx = {w:i for i,w in enumerate(Zl)}; n=len(Zl)
    rows=[]
    for u0,v0,u1,v1,_ in relations_2Bstar(d):
        if all(x in Z for x in (u0,v0,u1,v1)):
            r=[0]*n
            r[idx[u0]]+=1; r[idx[v0]]+=1; r[idx[u1]]-=1; r[idx[v1]]-=1
            if any(r): rows.append(r)
    if not rows: return n, []
    M=Matrix(rows)
    rank=M.rank()
    from sympy.matrices.normalforms import smith_normal_form
    S=smith_normal_form(M)
    tors=[S[i,i] for i in range(min(S.shape)) if S[i,i] not in (0,1,-1)]
    return n-rank, tors

for d in range(3,7):
    ws, good = enumerate_supports(d)
    print('d=%d' % d)
    seen={}
    for Z in good:
        if not Z: continue
        fr, tors = lattice_invariants(Z, d, ws)
        key=(len(Z), fr, tuple(tors))
        seen[key]=seen.get(key,0)+1
    for (sz,fr,tors),cnt in sorted(seen.items()):
        tag=''
        if sz==len(ws): tag=' (full support)'
        print('   |Z|=%3d : %d supports ; value torus dim %d ; torsion %s%s' % (sz,cnt,fr,list(tors) if tors else 'none',tag))
    sys.stdout.flush()
