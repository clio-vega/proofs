"""Test the commutator on Maya sets sampled DIRECTLY (lem:realise says every window
configuration is a partition), not on small partitions -- the two-bead sector needs a
window of size g=e+f and small partitions do not realise one."""
import sys, random
sys.path.insert(0,'.')
from independent_maya import L, inM, word, apply_R, commutator, alpha
from itertools import product

def maya_from_window(S_window, lo):
    """M = {n < -L} u {lo+i : i in S_window}; S_window a set of offsets >=0."""
    return frozenset(lo+i for i in S_window)

def scan(e, W, f, Wb, trials=4000, seed=0):
    """random Maya states over a window of width e+f+4 starting at 0."""
    rng = random.Random(seed)
    g = e+f; width = g+4
    nonzero = 0; twobead = 0
    for _ in range(trials):
        S = set(i for i in range(width) if rng.random() < 0.5)
        M = maya_from_window(S, 0)
        c = commutator({M:1}, e, W, f, Wb)
        if c:
            nonzero += 1
            for T in c:
                if len(set(T) ^ set(M)) >= 4: twobead += 1; break
    return nonzero, twobead

def exhaustive(e, W, f, Wb, width):
    n = 0
    for bits in product((0,1), repeat=width):
        M = maya_from_window({i for i in range(width) if bits[i]}, 0)
        if commutator({M:1}, e, W, f, Wb): n += 1
    return n, 2**width
