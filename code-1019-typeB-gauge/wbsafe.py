"""The analytically cancelled form of Motegi's W^B, and At^B built on it.

W^B := prod_k (1+u_k^{-1})^{-1} * F|_{gamma=-1}.  The prefactor has a POLE at
u_k = -1 and line 3c has a ZERO there; substituting u_k = -1 into the literal
expression gives zoo*0 = nan.  Cancelling first, per permutation sigma with
v_j = u_sigma(j):

  prod_{k=1}^n u_k/(u_k+1) * T_sigma^{(gamma=-1)}
    = (prod_k u_k) * [lines 1,2a,2b,3a,3b] * prod_{j>r}(v_j+1) / prod_{j=1}^n (v_j+1)
    = (prod_k u_k) * [lines] / prod_{j<=r}(v_j+1),

and the surviving 1/(v_j+1) cancels line 1's (v_j^2-1) = (v_j-1)(v_j+1):

    TYPE B HEAD FACTOR IS (v_j - 1)   where   TYPE C's IS (v_j^2 - 1).

Verified against the engine's WX(...,'B') at generic u: see v6_minusone.py V6.0.
"""
import sympy as sp
from itertools import permutations
from math import factorial
from q484_engine import EB, q


def TB(N, n, r, us, w, I):
    """T_sigma for type B, cancelled form: head factor (v-1), not (v^2-1)."""
    t = sp.Integer(1)
    for j in range(r):
        ij = I[j]
        t *= (1 - q) * (us[j] - 1)
        for k in range(ij + 1, N + 1):
            t *= (us[j] - q * w[k - 1])
        for k in range(1, ij):
            t *= (us[j] - w[k - 1])
    for j in range(r):
        for k in range(j + 1, r):
            t *= (q * us[j] - us[k]) / (us[j] - us[k]) * (us[j] * us[k] - 1)
    for j in range(r, n):
        for k in range(N):
            t *= (us[j] - q * w[k])
    for j in range(r, n):
        for k in range(r):
            t *= (q * us[j] - us[k]) / (us[j] - us[k]) * (us[j] * us[k] - 1)
    for j in range(r, n):
        for k in range(j + 1, n):
            t *= (us[j] * us[k] - q)
    return t


def WB_safe(N, n, r, u, w, I):
    tot = sp.Integer(0)
    for sg in permutations(range(n)):
        tot += TB(N, n, r, [u[sg[j]] for j in range(n)], w, I)
    pref = sp.Integer(1)
    for k in range(n):
        pref *= u[k] ** (-(N + n - 1))
    return sp.together(pref * tot / factorial(n - r))


def AtB(N, r, u, w, I):
    n = len(u)
    if n == 0:
        return sp.cancel(1 / EB(N, 0, 0, [], w))
    return sp.cancel(sp.together(WB_safe(N, n, r, u, w, I) / EB(N, n, r, u, w)))


def z(e):
    return sp.cancel(sp.together(e)) == 0
