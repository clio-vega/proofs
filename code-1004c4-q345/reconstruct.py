"""Attempt to RECONSTRUCT the 'increasing labeled chain' rule of
Lenart-Sottile math/0202090 Thm 2 from the constraint it must satisfy:

    sum over increasing labeled chains in [e,v]  of  x^{labels}  =  S_v .

I hold math/0202090 at extraction level 'abstract' only, so this is a derivation
attempt, not a reading.  Search space: (label interval) x (order condition).
"""
from itertools import permutations, product
from schub import length, schubert_table
from chains import chains

N = None

def interval(a, b, kind, n):
    if kind == "1..a":      return range(1, a+1)
    if kind == "1..b-1":    return range(1, b)
    if kind == "a..b-1":    return range(a, b)
    if kind == "a..a":      return range(a, a+1)
    if kind == "b-1..b-1":  return range(b-1, b)
    if kind == "1..n-1":    return range(1, n)
    raise ValueError(kind)

def key(a, b, kind):
    if kind == "a": return a
    if kind == "b": return b
    if kind == "len": return b-a
    if kind == "none": return 0
    raise ValueError(kind)

def ok_seq(qs, abs_, cond, kkind):
    if cond == "free":   return True
    if cond == "weak":   return all(qs[i] <= qs[i+1] for i in range(len(qs)-1))
    if cond == "strict": return all(qs[i] <  qs[i+1] for i in range(len(qs)-1))
    if cond == "weak-dec":   return all(qs[i] >= qs[i+1] for i in range(len(qs)-1))
    if cond == "strict-dec": return all(qs[i] >  qs[i+1] for i in range(len(qs)-1))
    if cond == "bjs":   # strict when key increases, weak otherwise  (BJS compatible seq)
        for i in range(len(qs)-1):
            ki, kj = key(*abs_[i], kkind), key(*abs_[i+1], kkind)
            if ki < kj:
                if not qs[i] < qs[i+1]: return False
            else:
                if not qs[i] <= qs[i+1]: return False
        return True
    if cond == "bjs-rev":
        for i in range(len(qs)-1):
            ki, kj = key(*abs_[i], kkind), key(*abs_[i+1], kkind)
            if ki > kj:
                if not qs[i] < qs[i+1]: return False
            else:
                if not qs[i] <= qs[i+1]: return False
        return True
    raise ValueError(cond)

def gen(v, n, ikind, cond, kkind, adj_only):
    """generating polynomial as dict: sorted exponent vector -> count"""
    e = tuple(range(1, n+1))
    out = {}
    for ch in chains(e, v, n):
        if adj_only and any(b != a+1 for (a, b) in ch):
            continue
        ranges = [interval(a, b, ikind, n) for (a, b) in ch]
        for qs in product(*ranges):
            if ok_seq(qs, ch, cond, kkind):
                ex = [0]*n
                for q in qs:
                    ex[q-1] += 1
                t = tuple(ex)
                out[t] = out.get(t, 0)+1
    return {k: c for k, c in out.items() if c}

def schub_dict(f, n):
    out = {}
    for ex, c in f.items():
        out[tuple(ex[:n])] = c
    return {k: c for k, c in out.items() if c}

if __name__ == "__main__":
    IK = ["1..a", "1..b-1", "a..b-1", "a..a", "b-1..b-1"]
    CD = ["free", "weak", "strict", "weak-dec", "strict-dec", "bjs", "bjs-rev"]
    KK = ["a", "b", "len"]
    survivors = []
    for n in (3,):
        tab = schubert_table(n, n)
        targets = {v: schub_dict(tab[v], n) for v in permutations(range(1, n+1))}
        for ikind in IK:
            for cond in CD:
                for kkind in (KK if cond.startswith("bjs") else ["none"]):
                    for adj in (False, True):
                        good = all(gen(v, n, ikind, cond, kkind, adj) == targets[v]
                                   for v in targets)
                        if good:
                            survivors.append((ikind, cond, kkind, adj))
    print("survivors at n=3:", survivors)
    # now test survivors at n=4 and n=5
    for (ikind, cond, kkind, adj) in survivors:
        for n in (4, 5):
            tab = schubert_table(n, n)
            targets = {v: schub_dict(tab[v], n) for v in permutations(range(1, n+1))}
            bad = [v for v in targets if gen(v, n, ikind, cond, kkind, adj) != targets[v]]
            print(f"  rule {(ikind,cond,kkind,adj)} at n={n}: "
                  f"{len(targets)-len(bad)}/{len(targets)} permutations correct"
                  + (f"  first failure {bad[0]}" if bad else ""))
