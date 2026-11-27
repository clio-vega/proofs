"""S1: is R^E_{k+3} (created by T_E's (0,b) at the exit R) of a packed shape for some nest?  Which classkey?
W1: tiles of E_{k+3} (each candidate type at class 1) vs G(L_0), G(L_1), G(E^L_j) j=2..k+2 over all corridor words."""
import itertools, sys
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
U = {'B': (-1, 1, 0, 0), 'C': (1, 0, 0, 0), 'A': (0, -1, 0, 0)}; W = {'B': (0, 0, 1, -1), 'C': (0, 0, -1, 0), 'A': (0, 0, 0, 1)}
add = lambda P, v: tuple(a + b for a, b in zip(P, v))
def packed_shape(R, T):
    S = set(R[1]); return any({p, add(p, U[T]), add(p, W[T]), add(add(p, U[T]), W[T])} == S for p in S)
aT = {'B': 'a1', 'C': 'a0'}
def pred(T, R, t):
    """the rhombus Rp with step(T,Rp,t) creating R (translation)"""
    rem, Rn, cre = step(T, R, t); d = tuple(a - b for a, b in zip(min(R[1]), min(Rn[1])))
    Rp = tr(R, d); assert step(T, Rp, t)[1] == R; return Rp
for TL in 'BC':
    TE = OTHER[TL]; R = canon(TE, 0)
    remEb, RE3, creEb = step(TE, R, 'b')
    print('S1 T_L', TL, 'R^E_{k+3} class', cls(TE, RE3), 'packed shape for', [T for T in 'ABC' if packed_shape(RE3, T)],
          '| classkey==K1[T_E]', classkey(RE3) == K1[TE], '| types at R^E_{k+3}:', sorted(configs_at(TE, RE3)))
    # also: L's created rhombus
    remLb, RL3, creLb = step(TL, R, 'b')
    print('    R^L_{k+3} packed shape for', [T for T in 'ABC' if packed_shape(RL3, T)], 'class', cls(TL, RL3))
    # and the canonical packed rhombus of T_E: its class / classkey, for comparison
