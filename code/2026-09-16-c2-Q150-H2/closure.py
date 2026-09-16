"""Theorem D: Z is admissible iff closed under the linked-pair flip rule.
   Hence the admissible sets form a closure system."""
import sys, itertools; sys.path.insert(0,'.')
from twobead import words, relations_2Bstar, support_ok
from supports import enumerate_supports

def linked_rule_closure(X, d):
    """smallest Z >= X closed under: u,v in Z delta-linked  =>  u^(delta), v^(n-delta+1) in Z"""
    n=d-1; Z=set(X)
    changed=True
    while changed:
        changed=False
        for u in list(Z):
            for v in list(Z):
                for delta in range(1,n+1):
                    if u[delta:]==v[:n-delta] and u[delta-1]==v[n-delta]:
                        a=list(u); a[delta-1]^=1; a=tuple(a)
                        b=list(v); b[n-delta]^=1; b=tuple(b)
                        if a not in Z: Z.add(a); changed=True
                        if b not in Z: Z.add(b); changed=True
    return frozenset(Z)

for d in (3,4,5,6):
    ws,good=enumerate_supports(d); rels=list(relations_2Bstar(d)); G=set(good)
    # (a) closed-under-rule  ==  admissible ?
    allsets = [frozenset(S) for r in range(len(ws)+1) for S in itertools.combinations(ws,r)] if d<=4 else None
    mism=0
    test = allsets if allsets is not None else list(G)+[frozenset(x) for x in itertools.combinations(ws,2)]
    for Z in test:
        closed = (linked_rule_closure(Z,d)==Z)
        adm = support_ok(Z,rels)
        if closed!=adm: mism+=1
    print('d=%d : rule-closed == admissible on %d test sets, mismatches %d' % (d,len(test),mism))
    # (b) intersection closure
    bad=sum(1 for a,b in itertools.combinations(good,2) if frozenset(a&b) not in G)
    print('      admissible sets closed under intersection: %s (%d pairs)' % (bad==0, bad))
    # (c) closure of a singleton
    from collections import Counter
    cnt=Counter()
    for w in ws:
        cl=linked_rule_closure({w},d)
        cnt[len(cl)]+=1
    print('      |closure of a single word|:', dict(sorted(cnt.items())), ' (total words %d)'%len(ws))
    sys.stdout.flush()
