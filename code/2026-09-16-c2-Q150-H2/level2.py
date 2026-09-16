import sys, itertools; sys.path.insert(0,'.')
from supports import enumerate_supports
from analyse import Cset
def mid(Z,d): return frozenset(w[1:d-2] for w in Z)
for d in (5,6,7):
    ws,good=enumerate_supports(d); prop=[Z for Z in good if 0<len(Z)<len(ws)]
    C=Cset(d,1,0); m=d-3
    Ss=[mid(Z,d) for Z in prop if Z<=C]
    full=frozenset(itertools.product((0,1),repeat=m))
    proper=[S for S in Ss if S!=full]
    maximal=[S for S in proper if not any(S<T for T in proper)]
    print('d=%d m=%d : %d admissible S, %d proper, %d MAXIMAL proper'%(d,m,len(Ss),len(proper),len(maximal)))
    print('    maximal sizes:', sorted(len(S) for S in maximal), ' (|full|=%d)'%len(full))
    # is each maximal a co-singleton?  a sub-cylinder {x_i=a, x_{m+1-i}=1-a}?
    cos=sum(1 for S in maximal if len(S)==len(full)-1)
    subcyl=0
    for S in maximal:
        for i in range(1,m+1):
            for a in (0,1):
                if 2*i!=m+1 and S=={x for x in full if x[i-1]==a and x[m-i]==1-a}: subcyl+=1
    print('    co-singletons: %d/%d ; sub-cylinders: %d/%d'%(cos,len(maximal),subcyl,len(maximal)))
    if m<=3:
        for S in maximal: print('      ', sorted(''.join(map(str,x)) for x in S))
    sys.stdout.flush()
