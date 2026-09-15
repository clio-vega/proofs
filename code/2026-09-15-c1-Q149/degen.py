"""Nonvanishing lemma, tested at e>=3: can sigma vanish somewhere without
one of the two operators being zero?  For each proper zero-pattern of W,
substitute and ask whether the ideal forces Wbar = 0."""
import sys, itertools
import sympy as sp
from solver import commutator_eqs_window, weight_syms
from engine import all_words

def test(e, f):
    uw = all_words(e-1); yw = all_words(f-1)
    Ws = weight_syms('W', e); Vs = weight_syms('V', f)
    eqs,_,_ = commutator_eqs_window(e, f, Lw=e+f+1)
    eqs = list(eqs)
    survivors = []
    for r in range(1, len(uw)):           # proper nonempty zero-sets
        for Z in itertools.combinations(uw, r):
            sub = {Ws[u]: 0 for u in Z}
            E = [sp.expand(q.subs(sub)) for q in eqs]
            E = [q for q in E if q != 0]
            vs = [Vs[y] for y in yw]
            # does the ideal force every Vbar to be nilpotent (i.e. Wbar == 0)?
            G = sp.groebner(E + [sp.Symbol('T')*sp.prod([Ws[u] for u in uw if u not in Z]) - 1],
                            *( [Ws[u] for u in uw if u not in Z] + vs + [sp.Symbol('T')] ),
                            order='lex')
            # Wbar == 0 forced  <=>  each V_y is in the radical
            forced = all(sp.groebner(list(G.exprs)+[sp.Symbol('S')*v-1],
                          *( [Ws[u] for u in uw if u not in Z] + vs
                             + [sp.Symbol('T'), sp.Symbol('S')]), order='lex').exprs == [sp.Integer(1)]
                         for v in vs)
            if not forced:
                survivors.append(Z)
    print(f'  (e,f)=({e},{f}): {2**(e-1)} words; proper zero-patterns with a '
          f'NONZERO Wbar solution: {len(survivors)}')
    for Z in survivors[:5]: print('     ', Z)
    return survivors

if __name__ == '__main__':
    test(3,4)
    test(3,5)
