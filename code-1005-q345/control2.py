"""Round 2 of planted controls, after two of three failed to fire.

Diagnosis of the two that didn't fire:
 (a) "non-strict lex" is VACUOUS: a label (k,b) cannot repeat in a chain.
     b=v(i) sits at position i<=k before the step and at position j>k after it,
     so for a second edge with the same b we would need i'=j>k, contradicting
     i'<=k.  So strict and non-strict lex give the same chains.  PROVED below
     by exhaustive search for a repeated label.
 (b) "b=u(j)" is the OTHER convention LS explicitly endorse (l.~356: "all the
     constructions and statements still hold if we uniformly set b=u(j)=w(i)").
     Its not firing CONFIRMS that remark; it is not a missed bug.

So those two ingredients need DIFFERENT controls.  Five real ones here.
"""
from collections import Counter
from itertools import permutations
from schub import length, rmul_t, schubert_table
import ls_side

ORIG_edges = ls_side.labelled_edges
ORIG_chains = ls_side.increasing_chains

def prop3(n):
    N = n+1
    w0 = tuple(range(n,1-1,-1))[:n]
    w0 = tuple(range(n,0,-1))
    delta = tuple(n-1-i for i in range(n-1))
    tab = schubert_table(n, N)
    bad = 0
    for u in permutations(range(1,n+1)):
        poly = Counter()
        for alpha,c in ls_side.I_table(u,w0,n).items():
            e = tuple([delta[i]-alpha[i] for i in range(n-1)]+[0]*(N-(n-1)))
            if min(e) < 0: return "FIRES (negative exponent)"
            poly[e]+=c
        if {e:c for e,c in poly.items() if c} != {e:c for e,c in tab[u].items() if c}:
            bad+=1
    return f"{bad} mismatches" + ("  <-- FIRES" if bad else "  <-- silent")

def chainfactory(cond):
    def f(u,w,n):
        target = length(w)
        def rec(v,last):
            if v==w: yield []; return
            if length(v)>=target: return
            for (k,b,v2) in ls_side.labelled_edges(v,n):
                if not cond((k,b),last): continue
                for rest in rec(v2,(k,b)): yield [(k,b)]+rest
        yield from rec(u,(0,0))
    return f

n = 4
# --- first: is a repeated label even possible? ---
rep = 0; tot = 0
for u in permutations(range(1,n+1)):
    for w in permutations(range(1,n+1)):
        if length(w)<=length(u): continue
        ls_side.increasing_chains = chainfactory(lambda a,b: True)   # ALL chains
        for ch in ls_side.increasing_chains(u,w,n):
            tot+=1
            if len(set(ch))!=len(ch): rep+=1
ls_side.increasing_chains = ORIG_chains
print(f"repeated-label check: {rep} of {tot} labelled chains (n=4) repeat a label  "
      f"-> non-strict lex is {'VACUOUSLY equal to strict' if rep==0 else 'a real weakening'}")

print("baseline                                        :", prop3(n))
tests = [
  ("no lex condition at all (count every chain)",      chainfactory(lambda a,b: True)),
  ("order by FIRST COORDINATE only, weakly increasing", chainfactory(lambda a,b: a[0]>=b[0])),
  ("lex on (b,k) instead of (k,b)",                     chainfactory(lambda a,b: (a[1],a[0])>(b[1],b[0]))),
  ("lex DEcreasing",                                    chainfactory(lambda a,b: b==(0,0) or a<b)),
  ("second coordinate must be weakly incr. in blocks",  chainfactory(lambda a,b: a[0]>b[0] or (a[0]==b[0] and a[1]>=b[1]))),
]
for desc,f in tests:
    ls_side.increasing_chains = f
    print(f"control: {desc:50s}:", prop3(n))
    ls_side.increasing_chains = ORIG_chains
print("baseline restored                               :", prop3(n))
