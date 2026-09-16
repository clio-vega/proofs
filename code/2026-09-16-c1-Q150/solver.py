"""
Commutator equations for shape-weighted ribbon operators, derived from the RAW
operator on Maya sets (no lemma input).

A Maya set is stored as a frozenset  S  of sites >= CUT, with the convention that
every site < CUT is a bead.  Initial states: S = [CUT, -L) u T with T an L-subset
of [-L, L).  Those are exactly the partitions inside an LxL box.  CUT is low
enough (-L-2g) that no bead below CUT has a legal move.
"""
import sympy as sp
from itertools import combinations
from engine import all_words, word_to_comp, comp_to_word

def sym(name, e, u):
    return sp.Symbol(f"{name}{e}_" + ("e" if not u else "".join(map(str,u))))

def weight_syms(name, e):
    return {u: sym(name, e, u) for u in all_words(e-1)}

def legal_moves(S, CUT, e):
    """all (b, word) for legal e-moves out of the Maya set (S, CUT)."""
    out = []
    lo = min(S) if S else CUT
    for b in sorted(S):
        if (b+e) in S: continue
        u = tuple(1 if (b+i) in S else 0 for i in range(1, e))
        out.append((b, u))
    return out

def apply_R(vec, e, W, CUT):
    out = {}
    for S, c in vec.items():
        for b, u in legal_moves(S, CUT, e):
            S2 = (S - {b}) | frozenset([b+e])
            out[S2] = out.get(S2, 0) + c*W[u]
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

def commutator_eqs(e, f, L=None, W=None, Wb=None):
    """Returns (set of distinct nonzero polynomial equations, n_states, n_elts)."""
    g = e+f
    if L is None: L = g+1
    CUT = -L-2*g
    if W is None:  W  = weight_syms("W", e)
    if Wb is None: Wb = weight_syms("V", f)
    base = frozenset(range(CUT, -L))
    eqs = set(); nstates = 0; nelts = 0
    for T in combinations(range(-L, L), L):
        S0 = base | frozenset(T); nstates += 1
        v = {S0: sp.Integer(1)}
        a = apply_R(apply_R(v, f, Wb, CUT), e, W, CUT)   # R_e R_f
        b = apply_R(apply_R(v, e, W,  CUT), f, Wb, CUT)  # R_f R_e
        for k in set(a) | set(b):
            val = sp.expand(a.get(k, 0) - b.get(k, 0))
            nelts += 1
            if val != 0:
                eqs.add(sp.Poly(val, *sorted(val.free_symbols, key=str)).as_expr()
                        if val.free_symbols else val)
    return eqs, nstates, nelts

def canon(eqs):
    """dedupe up to overall sign"""
    out = {}
    for q in eqs:
        p = sp.expand(q)
        key = min(sp.srepr(p), sp.srepr(sp.expand(-p)))
        out[key] = p
    return sorted(out.values(), key=lambda p: (len(p.free_symbols), sp.srepr(p)))


# ---------------- window enumeration (much cheaper than all charge-0 states) ----
def commutator_eqs_window(e, f, Lw=None, W=None, Wb=None):
    """States  M = Z_{<0} u T,  T subset of [0,Lw).  Every local window
    configuration of length <= Lw occurs; charge is irrelevant to the operator.
    Returns (eqs, n_states, n_elements)."""
    g = e+f
    if Lw is None: Lw = 2*f+2
    CUT = -2*g-2
    if W is None:  W  = weight_syms("W", e)
    if Wb is None: Wb = weight_syms("V", f)
    base = frozenset(range(CUT, 0))
    eqs = set(); nst = 0; nel = 0
    for mask in range(1 << Lw):
        T = frozenset(i for i in range(Lw) if (mask >> i) & 1)
        S0 = base | T; nst += 1
        v = {S0: sp.Integer(1)}
        a = apply_R(apply_R(v, f, Wb, CUT), e, W, CUT)
        b = apply_R(apply_R(v, e, W,  CUT), f, Wb, CUT)
        for k in set(a) | set(b):
            val = sp.expand(a.get(k, 0) - b.get(k, 0)); nel += 1
            if val != 0:
                eqs.add(sp.Poly(val, *sorted(val.free_symbols, key=str)).as_expr()
                        if val.free_symbols else val)
    return eqs, nst, nel
