"""Census for the gauge theorem.  Type C (h=2) unless stated.

Four blocks:
  A. the PURE adjunction lemma  Atil(R u D) = prod phi(y) Atil(R)  with the
     alphabet w FREE and (R,D) arbitrary -- this separates the lemma from the
     geometry, and is the only way d=4 is affordable.
  B. the GEOMETRIC instances: R,D drawn from w on the locus H_alpha, every root
     family, every M, every I.
  C. PLANTS, all drawn from OUTSIDE the hypothesis class (D not inverse-closed,
     or the prediction evaluated at the wrong h).
  D. the two algebraic identities the uniform proof turns on.
"""
import sympy as sp, time, itertools, sys
from gauge_check import *

t0 = time.time()
x, y, r1, r2 = sp.symbols('x y r1 r2')
score = {}
def rec(tag, res):
    if res is None: return
    if 'VACUOUS' in res:
        score.setdefault(tag + ' [vacuous]', [0, 0])
        s = score[tag + ' [vacuous]']; s[1] += 1; s[0] += int(res['VACUOUS']); return
    for k, v in res.items():
        if v is None: continue
        score.setdefault(f'{tag}/{k}', [0, 0])
        s = score[f'{tag}/{k}']; s[1] += 1; s[0] += int(bool(v))

# ============================ BLOCK A: pure adjunction lemma, free alphabet
print("--- A: pure adjunction lemma, alphabet free ---", flush=True)
EPS = [sp.Integer(1), sp.Integer(-1)]
DSETS = [("d=0 empty", []),
         ("d=1 fixed-pt eps=+1", [EPS[0]]),
         ("d=1 fixed-pt eps=-1", [EPS[1]]),
         ("d=2 free orbit {x,1/x}", [x, 1/x]),
         ("d=2 two fixed pts {1,1}", [EPS[0], EPS[0]]),
         ("d=2 two fixed pts {-1,-1}", [EPS[1], EPS[1]]),
         ("d=3 {x,1/x,-1}", [x, 1/x, EPS[1]]),
         ("d=4 {x,1/x,x,1/x}", [x, 1/x, x, 1/x])]
for N in (1, 2):
    w = wvars(N)
    for Rname, R in [("R=[]", []), ("R=[r1]", [r1]), ("R=[r1,r2]", [r1, r2])]:
        for I in [tuple(s) for k in range(N + 1) for s in itertools.combinations(range(1, N + 1), k)]:
            if len(R) < len(I):      # outside the lemma's hypothesis |R| >= r
                continue
            for dn, D in DSETS:
                if len(R) + len(D) > 5:
                    continue
                rec(f'A {dn}', adjunction(N, R, D, w, I, 'C', 2))
    print(f"  N={N} done [{time.time()-t0:.0f}s]", flush=True)

# ============================ BLOCK B: geometric instances
print("--- B: geometric instances on H_alpha ---", flush=True)
def subsets(N):
    return [tuple(s) for k in range(N + 1) for s in itertools.combinations(range(1, N + 1), k)]

for N in (2, 3):
    wsym = wvars(N)
    # --- boundary root 2 e_m : locus w_m = eps
    for m in range(1, N + 1):
        for eps in EPS:
            wl = list(wsym); wl[m - 1] = eps
            for M in subsets(N):
                if m not in M: continue
                R = [wl[j - 1] for j in M if j != m]
                for I in subsets(N):
                    if len(R) + 1 > 4: continue
                    rec(f'B boundary 2e_m N={N}', adjunction(N, R, [eps], wl, I, 'C', 2))
    # --- pair root e_a + e_b : locus w_a w_b = 1
    for a, b in itertools.combinations(range(1, N + 1), 2):
        wl = list(wsym); wl[a - 1] = x; wl[b - 1] = 1 / x
        for M in subsets(N):
            if a not in M or b not in M: continue
            R = [wl[j - 1] for j in M if j not in (a, b)]
            for I in subsets(N):
                if len(R) + 2 > 4: continue
                rec(f'B pair e_a+e_b N={N}', adjunction(N, R, [x, 1 / x], wl, I, 'C', 2))
    # --- NEW: composite point w_a = w_b = eps  (lies on {w_a^2=1} AND {w_a w_b=1})
    for a, b in itertools.combinations(range(1, N + 1), 2):
        for eps in EPS:
            wl = list(wsym); wl[a - 1] = eps; wl[b - 1] = eps
            for M in subsets(N):
                if a not in M or b not in M: continue
                R = [wl[j - 1] for j in M if j not in (a, b)]
                for I in subsets(N):
                    if len(R) + 2 > 4: continue
                    rec(f'B composite w_a=w_b=eps N={N}', adjunction(N, R, [eps, eps], wl, I, 'C', 2))
    print(f"  N={N} done [{time.time()-t0:.0f}s]", flush=True)

# --- A-type root e_i - e_j : the d=0 instance, stated as the congruence
print("--- B': A-type congruence (d=0) ---", flush=True)
for N in (2, 3):
    wsym = wvars(N)
    for i, j in itertools.combinations(range(1, N + 1), 2):
        wl = list(wsym); wl[j - 1] = wsym[i - 1]       # the locus w_j = w_i
        for M in subsets(N):
            if len(set(M) & {i, j}) != 1: continue
            sM = tuple(sorted(set(M) ^ {i, j}))
            for I in subsets(N):
                if len(M) > 3: continue
                L = AtilU(N, [wl[p - 1] for p in M], wl, I, 'C')
                Rt = AtilU(N, [wl[p - 1] for p in sM], wl, I, 'C')
                tag = f"B' A-type N={N}"
                score.setdefault(tag, [0, 0]); s = score[tag]
                s[1] += 1; s[0] += int(eq(L, Rt))

# ============================ BLOCK C: plants (outside the hypothesis class)
print("--- C: plants ---", flush=True)
PLANTS = [("P1 D=[x,y] not inverse-closed", [x, y], 2),
          ("P2 D=[2]  (2^2 != 1)",           [sp.Integer(2)], 2),
          ("P3 D=[x,1/x] scored at h=1",     [x, 1 / x], 1),
          ("P4 D=[x,x] not inverse-closed",  [x, x], 2)]
N = 2; w = wvars(N)
for pn, D, hh in PLANTS:
    for Rname, R in [("R=[]", []), ("R=[r1]", [r1])]:
        for I in [(), (1,)]:
            if len(R) < len(I): continue
            res = adjunction(N, R, D, w, I, 'C', hh, want_F=False, want_E=False)
            if res is None or 'VACUOUS' in res: continue
            score.setdefault(f'C {pn} MUST-FAIL', [0, 0])
            s = score[f'C {pn} MUST-FAIL']; s[1] += 1; s[0] += int(not res['Atil'])

# ============================ BLOCK D: the two algebraic identities
print("--- D: the identities the uniform proof turns on ---", flush=True)
def idents(D, must_hold):
    d = len(D)
    py = sp.Integer(1)
    for v in D: py *= v
    i1 = eq(P2(D), py**(d - 1) * Q2(D))
    den = sp.Integer(1)
    for v in D: den *= (1 - q * v**2)
    i2 = eq(Q2(D)**2 / Q1(D), 1 / den) if Q1(D) != 0 else None
    for nm, val in (('P2=(prod y)^{d-1} Q2', i1), ('Q2^2/Q1 = 1/prod(1-q y^2)', i2)):
        if val is None: continue
        tag = f'D {nm}' + ('' if must_hold else ' MUST-FAIL')
        score.setdefault(tag, [0, 0]); s = score[tag]
        s[1] += 1; s[0] += int(bool(val) == must_hold)
for _, D in DSETS:
    if D: idents(D, True)
idents([x, 1 / x, y, 1 / y], True)
idents([x, y], False)          # plant: not inverse-closed
idents([x, x], False)          # plant: not inverse-closed

# ============================ report
print(f"\n================ RESULTS  [{time.time()-t0:.0f}s] ================")
bad = 0
for k in sorted(score):
    g, t = score[k]
    flag = '' if g == t else '   <<<< DEFECT'
    if g != t: bad += 1
    print(f"  {g:4d}/{t:4d}  {k}{flag}")
print(f"\n{'ALL GREEN' if bad == 0 else str(bad) + ' DEFECTIVE LINES'}")
