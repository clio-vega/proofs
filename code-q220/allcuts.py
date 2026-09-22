"""Is the excess exactly the CUTS MS may not use?

For any x in Z/n define the window (x, x+n]: letters at position x sit before
the cut.  MS require x not in S u T (so the cut is clean); the pairing rule
itself makes sense for every cut.  Compute B_all = { b_x : x in Z/n }.
"""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *

def pairing_any_cut(S,T,x,n):
    """Bracket matching of the window (x, x+n]; ')' before '(' at a shared slot.
    Returns list of unmatched '(' positions (residues) in window order."""
    seq=[]
    for k in range(1,n+1):
        i=(x+k)%n
        if i in T: seq.append((')',i))
        if i in S: seq.append(('(',i))
    stack=[]; unmatched=[]
    for sym,i in seq:
        if sym=='(': stack.append(i)
        else:
            if stack: stack.pop()          # ')' closes the nearest open '('
    # recompute properly: '(' unmatched iff never closed
    stack=[]; closed=set()
    for sym,i in seq:
        if sym=='(': stack.append(i)
        else:
            if stack: closed.add(stack.pop())
    unmatched=[i for sym,i in seq if sym=='(' and i not in closed]
    return unmatched

def b_any_cut(S,T,x,n):
    un=pairing_any_cut(S,T,x,n)
    return un[0] if un else None

def B_all(S,T,n):
    return {b_any_cut(S,T,x,n) for x in range(n)} - {None}

# first: check pairing_any_cut reproduces ms_pairing when x not in S u T
for n in range(3,8):
    for S,T,_ in additive_pairs(n):
        for x in X_set(S,T,n):
            L1,_,_=ms_pairing(S,T,x,n)
            assert set(pairing_any_cut(S,T,x,n))==L1,(n,S,T,x)
print("cross-check: generic-cut bracket matching == MS pairing on legal cuts  [OK, n=3..7]")

tot=collections.Counter(); bad=[]
for n in range(3,8):
    eq=0; sub=0; N=0
    for S,T,_ in additive_pairs(n):
        E=clio_letters(S,T,n); Ba=B_all(S,T,n); N+=1
        if Ba<=E: sub+=1
        if Ba==E: eq+=1
        else: bad.append((n,S,T,E,Ba))
    print(f"n={n}: {N} additive pairs; B_all subset E: {sub}; B_all == E: {eq}")
    tot['N']+=N; tot['eq']+=eq; tot['sub']+=sub
print("TOTAL",dict(tot))
for b in bad[:8]: print("  MISMATCH n=%d S=%s T=%s E=%s B_all=%s"%(b[0],sorted(b[1]),sorted(b[2]),sorted(b[3]),sorted(b[4])))
