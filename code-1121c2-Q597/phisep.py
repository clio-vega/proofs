"""PROVE 11-21 c2: (Sh3-future) by a one-sided functional.  phi(P) = a0 + b0 (= eta_C of the Lambda_0 part), g(R) = phi(base R).
Facts read off the configuration table (Lemma word): along every B-step Delta g in {0,1}; along every C-step Delta g in {0,-1};
the class-0 b-step has Delta g = +1 for B and 0 for C.  For removed tiles t of a step from R: off(t) = phi(base t) - g(R).
T=B (A to B, J to C): t removed by A after the split has phi(base t) >= g(R*) + 1 + min off_B; removed by J after e has
  phi(base t) <= g(R*) + max off_C.  Contradiction iff 1 + min off_B > max off_C.
T=C (A to C, J to B), phi' = -phi: g' non-decreasing along C, J's split strict (-1 for phi'):
  contradiction iff min off'_C > max off'_B - 1.
Also prints the table, and the same for the other functionals in the 3-separator family (multiples).  usage: python3 phisep.py"""
import os, sys
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '.')
exec(open('words.py').read().split('\ndef run(T, L, st):')[0])
k0 = frozenset({(0, 1), (0, -1)})
def phi(P): return P[0] + P[1]
CF = {T: configs(T) for T in 'BC'}
dg, off = {}, {}
for T in 'BC':
    for (ck, typ), (Rc, rem, Rn, cre) in CF[T].items():
        dg[(T, ck, typ)] = phi(base(Rn)) - phi(base(Rc))
        for t in rem: off.setdefault(T, []).append((phi(base(t)) - phi(base(Rc)), ck, typ, t[0]))
for T in 'BC':
    print(T, 'Delta g by (class, type):', sorted((tuple(sorted(ck))[0], typ, str(v)) for (TT, ck, typ), v in dg.items() if TT == T))
    vals = [o[0] for o in off[T]]
    print(T, 'off(t) range over removed tiles:', str(min(vals)), str(max(vals)), sorted({(str(o[0]), o[3]) for o in off[T]}))
mB = min(o[0] for o in off['B']); MC = max(o[0] for o in off['C'])
print('T=B: 1 + min off_B =', 1 + mB, ' max off_C =', MC, ' -> separated:', 1 + mB > MC)
offpC = [-o[0] for o in off['C']]; offpB = [-o[0] for o in off['B']]
print("T=C: min off'_C =", min(offpC), " max off'_B - 1 =", max(offpB) - 1, ' -> separated:', min(offpC) > max(offpB) - 1)
