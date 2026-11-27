"""R1: M0 tiles removed by L_0 or L_1 (contact).  For each, either J_E removes it at E_0 or E_1, or its 4-offset
f = om(g, R) at J_E's rhombus R_2 is already PAST the removal box Q_{T_E} (psi-part > max, or A-part < min).
f is monotone along any terminating journey (Lemma window, G1c tex: psi-part non-decreasing, A-part non-increasing),
so J_E then never removes g at a step >= p_E+2.  Q_T is rebuilt from the local configuration table."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ETA = {'B': lambda a, b: b, 'C': lambda a, b: a + b}
def psi(P, T): return (min(-ETA[T](p[2], p[3]) for p in P), min(ETA[T](p[0], p[1]) for p in P))
def Af(P): return (min(p[0] for p in P), min(p[2] for p in P))
def om(t, R, T):
    a, b = psi(t[1], T), psi(R[1], T); c, d = Af(t[1]), Af(R[1]); return (a[0]-b[0], a[1]-b[1], c[0]-d[0], c[1]-d[1])
Q = {}
for T in 'BC':
    for c in (0, 1):
        R = canon(T, c)
        for typ, (rem, Rn, cre) in configs_at(T, R).items():
            for g in nonrh(rem): Q.setdefault(T, set()).add(om(g, R, T))
    print(T, '|Q_T| =', len(Q[T]))
def past(f, T):
    q = Q[T]
    return f[0] > max(x[0] for x in q) or f[1] > max(x[1] for x in q) or f[2] < min(x[2] for x in q) or f[3] < min(x[3] for x in q)
H = Counter()
for TL in 'BC':
    TE = OTHER[TL]
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, TE):
        cfL = configs_at(TL, RL1)
        if 'b' not in cfL: continue
        remL1, R2, creL1 = cfL['b']
        for tE1, (remE1, RE2, _) in configs_at(TE, RE1).items():
            if RE2 != R2: continue                 # E_1 must create R_2 (shadow: shared rhombus from step 2)
            GE01 = set(nonrh(HE)) | set(nonrh(remE1))
            for g in set(nonrh(HL)) | set(nonrh(remL1)):
                if g in GE01: H[(TL, tE1, 'removed by E_0/E_1')] += 1; continue
                f = om(g, R2, TE)
                H[(TL, tE1, 'past at R_2', past(f, TE), f if not past(f, TE) else '')] += 1
for k in sorted(H, key=str): print(H[k], k)
