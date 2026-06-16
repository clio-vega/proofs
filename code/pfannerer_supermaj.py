#!/usr/bin/env python3
"""
Job B (2026-06-16): Pfannerer super-maj vs the parity-twisted descent statistic s(T).

s(T) = sum_{i in Des(T)} w_i,  w_i = 2i-1 if (n-i) odd else 0.   G_lambda(q)=sum_T q^{s(T)}.

Pfannerer super-maj (per CODE.md, 2603.16598): signed tableau (Tp in SYT(lam), D subseteq {1..n}).
  super-descent i  iff  (i in Des(Tp) and i+1 not in D) or (i not in Des(Tp) and i in D).
  supermaj = sum of super-descent indices i.

Comparison 1: is there a sign-set rule D=D(Tp) making {supermaj} = {s(T)} as multisets?
Comparison 2: f^lam(q,t)=f^lam(q) prod_box (1+ t q^{c(box)}); set t=-1; do surviving
              contents match J* (= ord_{q=-1} G_lambda support / the d=4 box)?
"""
import sympy
from sympy import symbols, Poly
from itertools import combinations
from collections import Counter

# ---------- SYT generation + descents ----------
def syt_list(lam):
    target = list(lam)
    n = sum(lam)
    rowofval = [0] * (n + 1)   # rowofval[v] = row of value v, v=1..n
    results = []
    shape = [0] * (len(lam) + 1)
    def rec(v):
        if v > n:
            results.append(tuple(rowofval[1:]))  # rows of values 1..n, index 0..n-1
            return
        for r in range(len(target)):
            if shape[r] < target[r] and (r == 0 or shape[r] < shape[r - 1]):
                shape[r] += 1
                rowofval[v] = r
                rec(v + 1)
                shape[r] -= 1
    rec(1)
    return results, n

def descents(rowrec, n):
    # rowrec[k] = row of value k+1 (k=0..n-1). descent at i (1..n-1): value i+1 lower row than i.
    return set(i for i in range(1, n) if rowrec[i] > rowrec[i - 1])

def s_statistic(rowrec, n):
    s = 0
    for i in descents(rowrec, n):
        s += (2 * i - 1) if ((n - i) % 2 == 1) else 0
    return s

# ---------- Pfannerer super-maj ----------
def supermaj(Des, D, n):
    tot = 0
    for i in range(1, n):
        sd = (i in Des and (i + 1) not in D) or (i not in Des and i in D)
        if sd:
            tot += i
    return tot

def s_multiset(lam):
    sl, n = syt_list(lam)
    return Counter(s_statistic(rv, n) for rv in sl), n

# ---------- Comparison 1 ----------
def comparison1(lam):
    sl, n = syt_list(lam)
    smult = Counter(s_statistic(rv, n) for rv in sl)
    Deslist = [descents(rv, n) for rv in sl]

    # full super-maj distribution over ALL signed tableaux (Tp, D), D subseteq {1..n}
    full = Counter()
    for Des in Deslist:
        for k in range(0, n + 1):
            for Dt in combinations(range(1, n + 1), k):
                full[supermaj(Des, set(Dt), n)] += 1

    # candidate D-rules (one D per Tp):
    def rule_descent(Des):       # D = Des(Tp)
        return set(Des)
    def rule_parityflip(Des):    # D = {i : n-i odd}  (the cells where n-i parity is odd)
        return set(i for i in range(1, n + 1) if (n - i) % 2 == 1)
    def rule_empty(Des):
        return set()
    def rule_all(Des):
        return set(range(1, n + 1))
    def rule_nondescent_parity(Des):  # D = non-descents with n-i odd
        return set(i for i in range(1, n) if i not in Des and (n - i) % 2 == 1)

    rules = {'D=Des': rule_descent, 'D=parityflip': rule_parityflip,
             'D=empty': rule_empty, 'D=all': rule_all,
             'D=nondesc&parity': rule_nondescent_parity}
    rule_match = {}
    for name, rule in rules.items():
        mult = Counter(supermaj(Des, rule(Des), n) for Des in Deslist)
        rule_match[name] = (mult == smult, mult)

    # existence: for each Tp, can SOME D give value = some s-value, matched as multiset?
    # (greedy bipartite feasibility: does the set of achievable supermaj per Tp cover smult?)
    achievable = []
    for Des in Deslist:
        vals = set()
        for k in range(0, n + 1):
            for Dt in combinations(range(1, n + 1), k):
                vals.add(supermaj(Des, set(Dt), n))
        achievable.append(vals)
    exists = hall_match(achievable, smult)
    return smult, full, rule_match, exists, n

def hall_match(achievable, target_mult):
    """Can we assign each Tp a value from its achievable set so that the multiset = target?"""
    # expand target into list, do simple augmenting-path bipartite matching
    targets = []
    for v, c in target_mult.items():
        targets += [v] * c
    # build adjacency: tableau index -> list of target-slot indices it can fill
    adj = [[j for j, t in enumerate(targets) if t in ach] for ach in achievable]
    matchT = [-1] * len(targets)
    def try_assign(u, seen):
        for j in adj[u]:
            if not seen[j]:
                seen[j] = True
                if matchT[j] == -1 or try_assign(matchT[j], seen):
                    matchT[j] = u
                    return True
        return False
    cnt = 0
    for u in range(len(achievable)):
        if try_assign(u, [False] * len(targets)):
            cnt += 1
    return cnt == len(achievable) == len(targets)

# ---------- Comparison 2 ----------
def contents(lam):
    return [j - i for i, p in enumerate(lam) for j in range(p)]

def Gpoly(lam):
    sl, n = syt_list(lam)
    q = symbols('q')
    return Poly(sum(q ** s_statistic(rv, n) for rv in sl), q), q, n

def root_mult(poly, q, root):
    mult = 0; p = poly
    while not p.is_zero and p.eval(root) == 0 and mult < 500:
        quo, rem = sympy.div(p.as_expr(), (q - root), q)
        p = Poly(quo, q); mult += 1
    return mult

def comparison2(lam):
    q = symbols('q'); t = symbols('t')
    poly, q, n = Gpoly(lam)
    ordm1 = root_mult(poly, q, -1)
    # f^lam(q,t) = f^lam(q) * prod (1 + t q^{c}); at t=-1: prod (1 - q^{c})
    cs = contents(lam)
    prod = sympy.prod([(1 - q ** c) for c in cs])
    prod = sympy.expand(prod)
    # surviving contents under t=-1 cancellation: a factor (1-q^c) is "trivial" when c=0
    # (gives 1-1=0 only at q=1; but as polynomial 1-q^0 = 0 identically!).  Track c=0 cells.
    zero_content_cells = [c for c in cs if c == 0]
    # multiplicity of (q+1) i.e. q=-1 root in prod(1-q^c):
    Pp = Poly(prod, q) if prod != 0 else None
    ordm1_prod = root_mult(Pp, q, -1) if (Pp is not None and not Pp.is_zero) else None
    # content multiset
    cmult = Counter(cs)
    return {'lam': lam, 'n': n, 'G_at_-1_mult': ordm1, 'contents': sorted(cs),
            'content_mult': dict(sorted(cmult.items())),
            'zero_content_count': len(zero_content_cells),
            'prod(1-q^c)_is_zero': prod == 0,
            'prod_at_-1_mult': ordm1_prod, 'G_poly': poly.as_expr()}

# ---------- main ----------
if __name__ == "__main__":
    shapes = [(3, 1), (2, 1, 1), (2, 2), (3, 2), (2, 2, 1), (4, 2), (3, 1, 1), (3, 2, 1)]
    print("=" * 74)
    print("COMPARISON 1: Pfannerer super-maj vs s(T) multiset")
    print("=" * 74)
    for lam in shapes:
        smult, full, rule_match, exists, n = comparison1(lam)
        print(f"\nlambda={lam}  n={n}  |SYT|={sum(smult.values())}")
        print(f"  s(T) multiset:        {dict(sorted(smult.items()))}")
        anymatch = [name for name, (ok, _) in rule_match.items() if ok]
        for name, (ok, mult) in rule_match.items():
            flag = "  <== MATCH" if ok else ""
            print(f"  supermaj[{name}]: {dict(sorted(mult.items()))}{flag}")
        print(f"  EXISTS D-rule with multiset = s(T)?  {exists}")
        print(f"  (full supermaj dist over all signed tableaux, smallest 8: "
              f"{dict(sorted(full.items())[:8])})")

    print("\n" + "=" * 74)
    print("COMPARISON 2: t=-1 content cancellation vs ord_{q=-1} G_lambda")
    print("=" * 74)
    for lam in [(3, 1), (2, 1, 1), (2, 2), (3, 2), (3, 2, 1)]:
        r = comparison2(lam)
        print(f"\nlambda={lam}  n={r['n']}")
        print(f"  contents (multiset): {r['content_mult']}")
        print(f"  ord_(q=-1) G_lambda = {r['G_at_-1_mult']}")
        print(f"  #cells with content 0: {r['zero_content_count']}")
        print(f"  prod(1-q^c) identically 0? {r['prod(1-q^c)_is_zero']}  "
              f"ord_(q=-1) of prod = {r['prod_at_-1_mult']}")
        print(f"  G_lambda(q) = {r['G_poly']}")
