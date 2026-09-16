import sys, itertools; sys.path.insert(0,'.')
from twobead import words
from supports import enumerate_supports
from analyse import Cset

def middle(Z, d):
    return frozenset(w[1:d-2] for w in Z)

for d in (5,6,7):
    ws, good = enumerate_supports(d)
    prop=[Z for Z in good if 0<len(Z)<len(ws)]
    C=Cset(d,1,0)
    Ss=sorted((middle(Z,d) for Z in prop if Z <= C), key=lambda S:(len(S),sorted(S)))
    m=d-3
    print('d=%d  m=%d : %d admissible S of the %d nonempty subsets of {0,1}^%d'%(d,m,len(Ss),2**(2**m)-1,m))
    if m<=3:
        for S in Ss: print('     ', sorted(''.join(map(str,x)) for x in S))
    Sset=set(Ss)
    inter=all((not (a&c)) or frozenset(a&c) in Sset for a,c in itertools.combinations(Ss,2))
    union=all(frozenset(a|c) in Sset for a,c in itertools.combinations(Ss,2))
    comp=all(frozenset(tuple(1-t for t in x) for x in S) in Sset for S in Ss)
    rev =all(frozenset(tuple(reversed(x)) for x in S) in Sset for S in Ss)
    print('     closed under: intersection %s | union %s | letterwise-complement %s | reversal %s'%(inter,union,comp,rev))
    sys.stdout.flush()
