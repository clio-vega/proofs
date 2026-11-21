"""PROVE 11-21 c2: (Sh1) from the contact configuration (CC), finite derivation on the configuration table (Lemma word).
For each observed clean contact configuration (T, J types at c0, c0+1, A's type at p, offset lexmin R_p - lexmin R^o_c0):
 place J's two steps and A's step p;  check H_p & H^o_c0 = {sq}, H_p & H^o_c0+1 = {tri};
 the tiles on A's board near R_{p+1} include P = (non-rh tiles of H^o_c0, H^o_c0+1) minus H_p, and A's created tiles at p;
 for each class-1 type t of T at R_{p+1}: are all non-rh tiles of its good configuration in that set?  (expect: only b)
 then: image of A's b at R_{p+1} == R^o_{c0+2}, of class 0;  the non-rh tiles of A's b hexagon are J tiles of c0/c0+1.
Lesion: shift the offset by each unit vector (12): the derivation must fail (b-hexagon not available or image differs).
usage: python3 hsderive.py"""
import os, sys
from fractions import Fraction as Fr
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '.')
exec(open('words.py').read().split('\ndef run(T, L, st):')[0])
CF = {T: configs(T) for T in 'BC'}
k0 = frozenset({(0, 1), (0, -1)})
K1 = {'B': frozenset({(1, -1), (-1, 1)}), 'C': frozenset({(1, 0), (-1, 0)})}
OTHER = {'B': 'C', 'C': 'B'}
def step(T, R, typ):
    Rc, rem, Rn, cre = CF[T][(classkey(R), typ)]
    sh = base(R); assert trs(Rc, sh) == R
    return frozenset(trs(t, sh) for t in rem), trs(Rn, sh), frozenset(trs(t, sh) for t in cre)
def canon(T, cls):
    return CF[T][((k0 if cls == 0 else K1[T]), 'a0')][0]
CCS = [('B', (0, 'a0'), (0, 'a1'), (1, 'a1'), (0, 0, 0, 0)), ('B', (1, 'b'), (0, 'a1'), (1, 'a1'), (0, 0, 1, -1)),
       ('C', (0, 'a1'), (0, 'a0'), (1, 'a0'), (0, 1, 0, 1)), ('C', (1, 'b'), (0, 'a0'), (1, 'a0'), (1, 0, 0, 0))]
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
W12 = [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
def derive(T, t0, t1, tA, off, verbose=True):
    Tp = OTHER[T]
    RJ0 = canon(Tp, t0[0])
    H0, RJ1, C0 = step(Tp, RJ0, t0[1]); assert classkey(RJ1) == (k0 if t1[0] == 0 else K1[Tp])
    H1, RJ2, C1 = step(Tp, RJ1, t1[1])
    RAc = canon(T, tA[0]); sh = tuple(b + o - a for a, b, o in zip(base(RAc), base(RJ0), off)); RAp = trs(RAc, sh)
    Hp, RAp1, Cp = step(T, RAp, tA[1])
    nr = lambda S: {u for u in S if u[0] != 'rh'}
    kinds = lambda S: tuple(sorted(u[0] for u in S))
    r = {'Hp&H0': kinds(nr(Hp) & nr(H0)), 'Hp&H1': kinds(nr(Hp) & nr(H1))}
    avail = (nr(H0) | nr(H1)) - nr(Hp) | nr(Cp)
    ok_types = []
    for typ in ('a0', 'a1', 'b'):
        Ht, Rt, Ct = step(T, RAp1, typ)
        if nr(Ht) <= avail: ok_types.append(typ)
    r['class-1 types with hexagon present'] = tuple(ok_types)
    Hb, Rb, Cb = step(T, RAp1, 'b')
    r['b image == R^o_c0+2'] = Rb == RJ2
    r['b image class 0'] = classkey(Rb) == k0
    r['b hexagon tiles from J c0/c0+1 only'] = nr(Hb) <= (nr(H0) | nr(H1))
    r['OK'] = (r['Hp&H0'] == ('sq',) and r['Hp&H1'] == ('tri',) and 'b' in ok_types and r['b image == R^o_c0+2'] and r['b image class 0'] and r['b hexagon tiles from J c0/c0+1 only'])
    return r
for cc in CCS:
    r = derive(*cc); print(cc, r)
print('LESION: offset shifted by each of 12 unit vectors')
for cc in CCS:
    T, t0, t1, tA, off = cc; n_ok = 0
    for v in W12:
        try:
            r = derive(T, t0, t1, tA, tuple(a + b for a, b in zip(off, v)))
            n_ok += r['OK']
        except (KeyError, AssertionError):
            pass
    print(cc[:4], 'shifted offsets still deriving the handshake:', n_ok, 'of 12')
