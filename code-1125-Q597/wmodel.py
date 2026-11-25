"""shared word-model helpers (from wordlast.py)"""

import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597')
from collections import Counter
import words as W

CL0 = frozenset({(0, 1), (0, -1)})
CL1 = {'B': frozenset({(1, -1), (-1, 1)}), 'C': frozenset({(1, 0), (-1, 0)})}
AT = {'B': 'a1', 'C': 'a0'}; OTHER = {'B': 'C', 'C': 'B'}
CF = {T: W.configs(T) for T in 'BC'}
def vsub(P, Q): return tuple(a - b for a, b in zip(P, Q))
def vadd(P, Q): return tuple(a + b for a, b in zip(P, Q))
def step(T, R, typ):
    """(removed non-rh tiles, created non-rh tiles, image rhombus) of the good T-configuration of type typ at R"""
    Rc, rem, Rn, cre = CF[T][(W.classkey(R), typ)]
    sh = W.base(R); assert W.trs(Rc, sh) == R
    return (frozenset(W.trs(t, sh) for t in rem), frozenset(W.trs(t, sh) for t in cre), W.trs(Rn, sh))
def run_word(T, R0, word):
    """list of (R_k, removed_k, created_k) and final rhombus"""
    out = []; R = R0
    for (cl, typ) in word:
        ck = W.classkey(R); want = CL0 if cl == 0 else CL1[T]
        if ck != want: return None, None
        rem, cre, Rn = step(T, R, typ); out.append((R, rem, cre)); R = Rn
    return out, R
def canon(T):
    """canonical class-1 rhombus of T at base 0"""
    for (ck, typ), v in CF[T].items():
        if ck == CL1[T]: return v[0]
