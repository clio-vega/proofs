"""Q227: the bead (TASEP) model for Theorem 4.1's exchange move.

A bead configuration is A subset Z/N, 1 <= |A| <= N-1 (beads = particles on a
ring of N sites).  Lam's nilCoxeter action (math/0501335 sec.7), in the bead
coordinates derived in 2026-09-21-c2 (cylact.py docstring):

    u_i . A  =  (A \ {i}) u {i+1}     if i in A and i+1 not in A
             =  0                      otherwise.

This is exactly the TASEP exclusion rule: the particle at site i hops one step
clockwise iff site i+1 is empty.

Everything below is at the level of residues mod N; the rank of the cylindric
shape plays no role in whether an action is zero.
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-q220')
from affine import (u_S, window, compose, length, runs, clio_letters,
                    clio_moves, additive_pairs, ms_etilde, X_set)

N_STR = "bead model, residues only"


def u_act(A, i, n):
    """u_i . A.  A a frozenset of residues.  Returns None for 0."""
    i %= n
    if i not in A or (i + 1) % n in A:
        return None
    return frozenset((A - {i}) | {(i + 1) % n})


def word_cd(S, n):
    """The cyclically decreasing reduced word for u_S, LEFT to RIGHT:
    run [m,M] contributes s_M s_{M-1} ... s_m; runs concatenated."""
    w = []
    for r in runs(S, n):
        w.extend(reversed(r))          # r = [m, m+1, ..., M] -> M,...,m
    return w


def act_word(A, word, n, rightmost_first=True):
    """Apply the word to A.  rightmost_first=True means the RIGHTMOST letter
    of the word acts first (the standard 'leftmost acts last' convention of
    cylact.py)."""
    cur = A
    seq = reversed(word) if rightmost_first else word
    for i in seq:
        cur = u_act(cur, i, n)
        if cur is None:
            return None
    return cur


def act_cd(A, S, n, rightmost_first=True):
    return act_word(A, word_cd(S, n), n, rightmost_first)
