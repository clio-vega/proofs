#!/usr/bin/env python3
"""Reproduces every number quoted in section 7 of
proofs/2026-09-21-c1-affine-stanley-exchange.tex.   Run:  python3 verify_all.py
Only affstan.py (the direct implementation of (def-affine-stanley)) is imported."""
from itertools import permutations, combinations
from collections import Counter
from affstan import (identity, rmul_s, length, cyc_dec_elements,
                     affine_stanley_support, elements_of_length)

# ----------------------------------------------------------------- helpers
def ext(w, j, n):
    k, rem = divmod(j - 1, n)
    return w[rem] + k * n

def runs_Z(S, n):
    """maximal runs of the lift of S, one representative per period: (m, M), m in [0,n-1]"""
    out = []
    for m in range(n):
        if m in S and (m - 1) % n not in S:
            M = m
            while (M + 1) % n in S: M += 1
            out.append((m, M))
    return out

def reduced_words(w, n):
    if length(w, n) == 0: return 1
    return sum(reduced_words(rmul_s(w, i, n), n)
               for i in range(n) if length(rmul_s(w, i, n), n) == length(w, n) - 1)

def rmul_word(w, word, n):
    v = w
    for i in word:
        v2 = rmul_s(v, i, n)
        if length(v2, n) != length(v, n) + 1: return None, False
        v = v2
    return v, True

def two_factorisations(v, n):
    cd = cyc_dec_elements(n); lv = length(v, n); out = []
    for Sf, (uS, size, word) in cd.items():
        if size > lv: continue
        for Tf, (uT, size2, word2) in cd.items():
            if size + size2 != lv: continue
            z, ok = rmul_word(uS, word2, n)
            if ok and z == v: out.append((Sf, Tf))
    return out

def F(v, n): return sorted({len(S) for S, _ in two_factorisations(v, n)})

def is_M_convex(supp):
    S = set(supp); r = len(next(iter(S)))
    for a in S:
        for b in S:
            for i in range(r):
                if a[i] <= b[i]: continue
                if not any(a[j] < b[j] and
                           tuple(a[k] - (k == i) + (k == j) for k in range(r)) in S
                           for j in range(r)):
                    return False
    return True

R = {}
# --------------------------------------------------- anchor 1: reduced words
ok = bad = 0
for n in (2, 3, 4):
    byl = elements_of_length(n, 5)
    for L in range(1, 5):
        for w in byl[L]:
            c = affine_stanley_support(w, n, L).get(tuple([1] * L), 0)
            if c == reduced_words(w, n): ok += 1
            else: bad += 1
R['anchor: coeff(x_1..x_l) = #reduced words'] = (ok, bad)

# --------------------------------------------------- anchor 2: symmetry
ok = bad = 0
for n in (3, 4):
    byl = elements_of_length(n, 5)
    for L in range(1, 6):
        for w in byl[L]:
            for r in (2, 3, 4):
                Fw = affine_stanley_support(w, n, r)
                if all(Fw.get(p, 0) == c for a, c in Fw.items() for p in permutations(a)):
                    ok += 1
                else: bad += 1
R['anchor: F~_w symmetric'] = (ok, bad)

# --------------------------------------------------- block model
bad = 0; tot = 0
for n in (2, 3, 4, 5, 6):
    for Sf, (w, size, word) in cyc_dec_elements(n).items():
        G = {i for i in range(n) if i not in Sf}
        if not G: continue
        win = []
        for j in range(1, n + 1):
            if (j - 1) % n not in G: win.append(j - 1)
            else:
                c = j
                while c % n not in G: c += 1
                win.append(c)
        tot += 1
        if tuple(win) != w or n - len(G) != size: bad += 1
R['block model eq.(2) = reduced-word construction'] = (tot, bad)

# --------------------------------------------------- Lemma "add"
bad = 0; tot = 0
for n in (3, 4, 5):
    cd = cyc_dec_elements(n); byl = elements_of_length(n, 4)
    for L in range(5):
        for w in byl[L]:
            for Sf, (uS, size, word) in cd.items():
                if not Sf: continue
                z = w
                for i in word: z = rmul_s(z, i, n)
                additive = (length(z, n) == length(w, n) + size)
                crit = all(ext(w, M + 1, n) > max(ext(w, j, n) for j in range(m, M + 1))
                           for (m, M) in runs_Z(Sf, n))
                tot += 1
                if additive != crit: bad += 1
R['Lemma 2.3 (length-additivity criterion)'] = (tot, bad)

# --------------------------------------------------- Thm localisation
def fibre_pq(p, q, n):
    cd = cyc_dec_elements(n); lq = length(q, n); lp = length(p, n); vals = set()
    for Sf, (uS, size, word) in cd.items():
        x, ok1 = rmul_word(p, word, n)
        if not ok1: continue
        for Tf, (uT, s2, word2) in cd.items():
            if size + s2 != lq - lp: continue
            y, ok2 = rmul_word(x, word2, n)
            if ok2 and y == q: vals.add(size); break
    return sorted(vals)

bad = 0; tot = 0
for n in (3, 4):
    byl = elements_of_length(n, 5)
    for L in range(6):
        for p in byl[L]:
            for LL in range(6 - L):
                for vv in byl[LL]:
                    wd = []; u = vv
                    while length(u, n) > 0:
                        for i in range(n):
                            if length(rmul_s(u, i, n), n) == length(u, n) - 1:
                                wd.append(i); u = rmul_s(u, i, n); break
                    q, ok1 = rmul_word(p, wd[::-1], n)
                    if not ok1: continue
                    tot += 1
                    if fibre_pq(p, q, n) != F(vv, n): bad += 1
R['Thm 3.1 (fibre depends only on v = p^{-1}q)'] = (tot, bad)

# --------------------------------------------------- Thm exchange
def move(S, T, n):
    res = []
    for (m, M) in runs_Z(S, n):
        if m % n in T: continue
        e = next((j for j in range(m, M + 1) if (j + 1) % n not in T), None)
        if e is None: continue
        res.append((frozenset(S - {e % n}), frozenset(T | {m % n})))
    return res

tot = ok = invalid = missed = 0
step1_tot = step1_bad = 0
appl_tot = appl_bad = 0
id_tot = id_bad = 0
for n in (3, 4, 5, 6, 7):
    Lm = 5 if n < 7 else 4
    byl = elements_of_length(n, Lm); cd = cyc_dec_elements(n)
    for L in range(1, Lm + 1):
        for v in byl[L]:
            tf = two_factorisations(v, n)
            if not tf: continue
            allp = {(frozenset(S), frozenset(T)) for S, T in tf}
            for (S, T) in tf:
                for (m, M) in runs_Z(S, n):          # Lemma 4.2 Step 1
                    if m % n in T: continue
                    step1_tot += 1
                    if not all(j % n in T for j in range(m, M + 2)): step1_bad += 1
                if len(S) <= len(T): continue
                tot += 1
                appl_tot += 1
                if not [(m, M) for (m, M) in runs_Z(S, n)
                        if m % n not in T and any((j + 1) % n not in T for j in range(m, M + 1))]:
                    appl_bad += 1
                mv = move(S, T, n)
                for p in mv:
                    id_tot += 1
                    if p not in allp: id_bad += 1
                if any(p in allp for p in mv): ok += 1
                else: missed += 1
                invalid += sum(1 for p in mv if p not in allp)
# NB: step1 above is stated for m IN T; recompute correctly
step1_tot = step1_bad = 0
for n in (3, 4, 5, 6):
    byl = elements_of_length(n, 5)
    for L in range(1, 6):
        for v in byl[L]:
            for S, T in two_factorisations(v, n):
                for (m, M) in runs_Z(S, n):
                    if m % n not in T: continue
                    step1_tot += 1
                    if not all(j % n in T for j in range(m, M + 2)): step1_bad += 1
R['Thm 4.1 exchange rule yields a valid smaller factorisation'] = (tot, missed)
R['Thm 4.1 exchange rule never yields an invalid pair (= Lemma 4.3)'] = (id_tot, id_bad)
R['Lemma 4.2 Step 1 (m in T => [m,M+1] subset T)'] = (step1_tot, step1_bad)
R['Lemma 4.2 conclusion (|S|>|T| => usable run exists)'] = (appl_tot, appl_bad)

# --------------------------------------------------- (H4) and M-convexity
tot = noPM = noM = 0; widths = Counter(); suppsizes = Counter()
for n, Lmax in ((2, 7), (3, 6), (4, 5), (5, 4)):
    byl = elements_of_length(n, Lmax)
    for L in range(Lmax + 1):
        for w in byl[L]:
            for r in (2, 3, 4, 5):
                supp = list(affine_stanley_support(w, n, r).keys())
                if not supp: continue
                tot += 1; suppsizes[len(supp)] += 1
                best = [max(sum(a[:t]) for a in supp) for t in range(1, r + 1)]
                if not any(all(sum(a[:t]) == best[t - 1] for t in range(1, r + 1)) for a in supp):
                    noPM += 1
                if not is_M_convex(supp): noM += 1
R['(H4) simultaneous prefix maximum exists'] = (tot, noPM)
R['M-convexity of every support (verifies the code vs WZZ Thm 4.7)'] = (tot, noM)

for n in (4, 5):
    Lmax = 5 if n == 4 else 4
    byl = elements_of_length(n, Lmax); wd = Counter(); seen = set()
    for L in range(Lmax + 1):
        for v in byl[L]:
            if v in seen: continue
            seen.add(v)
            f = F(v, n)
            if f: wd[len(f)] += 1
    print(f"  fibre-width histogram n={n}: {dict(sorted(wd.items()))}")

print()
for k, (t, b) in R.items():
    print(f"  {'OK ' if b == 0 else 'FAIL'}  {k}: {t} tests, {b} failures")
print(f"\n  max support size seen: {max(suppsizes)}")

# --------------------------------------------------- Cor 4.4: the anti-automorphism Psi
def any_reduced_word(w, n):
    wd = []; u = w
    while length(u, n) > 0:
        for i in range(n):
            if length(rmul_s(u, i, n), n) == length(u, n) - 1:
                wd.append(i); u = rmul_s(u, i, n); break
    return wd[::-1]

def Psi(w, n):
    """phi(w^{-1}) with phi: s_i -> s_{-i}: the element of the negated, reversed word."""
    z = identity(n)
    for i in reversed(any_reduced_word(w, n)):
        z = rmul_s(z, (-i) % n, n)
    return z

bad_u = bad_f = tot_f = 0
for n in (3, 4, 5, 6):
    cd = cyc_dec_elements(n)
    for Sf, (uS, size, word) in cd.items():
        if not Sf: continue
        if Psi(uS, n) != cd[frozenset((-i) % n for i in Sf)][0]: bad_u += 1
    Lm = 5 if n < 6 else 4
    byl = elements_of_length(n, Lm)
    for L in range(1, Lm + 1):
        for v in byl[L]:
            f = F(v, n)
            if not f: continue
            tot_f += 1
            if sorted(L - k for k in f) != F(Psi(v, n), n): bad_f += 1
print(f"  {'OK ' if bad_u == 0 else 'FAIL'}  Cor 4.4: Psi(u_S) = u_(-S): 0 failures" if bad_u == 0
      else f"  FAIL  Cor 4.4: Psi(u_S) = u_(-S): {bad_u} failures")
print(f"  {'OK ' if bad_f == 0 else 'FAIL'}  Cor 4.4: F(Psi(v)) = l(v) - F(v): {tot_f} tests, {bad_f} failures")
