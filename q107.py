"""
Q107 sweep.  c(R,mu) := n(lambda) - n(mu) where lambda/mu is a border strip of shape R.

ENGINE 1 (here): abacus / Maya set.  Border strips are bead moves p -> p+e.
Predicted (hand-derived this session):
      c  =  (r-1)*e + maj(R)
where r = topmost row of the strip (1-indexed) and maj(R) = sum of the descent
set of the strip read from its BOTTOM-LEFT end.
"""
from itertools import combinations

def n_stat(lam):
    return sum(i*p for i, p in enumerate(lam))     # sum (i-1)lam_i, 0-indexed i

def maya(lam, L):
    """Maya set truncated: beads {lam_j - j : j=1..L}, shifted by L so >=0."""
    lam = list(lam) + [0]*(L-len(lam))
    return frozenset(lam[j] - (j+1) + L for j in range(L))

def from_maya(M, L):
    bs = sorted(M, reverse=True)
    lam = [bs[j] - L + (j+1) for j in range(len(bs))]
    return tuple(x for x in lam if x > 0)

def parts_upto(n):
    out = []
    def rec(rem, mx, acc):
        out.append(tuple(acc))
        for k in range(min(rem, mx), 0, -1):
            rec(rem-k, k, acc+[k])
    rec(n, n, [])
    return sorted(set(out))

def strip_data(mu, lam):
    """rows occupied, row-lengths top->bottom, top row r (1-indexed)."""
    mul = list(mu) + [0]*(len(lam)-len(mu))
    rows = [(i+1, lam[i]-mul[i]) for i in range(len(lam)) if lam[i] > mul[i]]
    r = rows[0][0]
    a = [x for _, x in rows]                      # top -> bottom
    return r, tuple(a)

def maj_of_strip(a):
    """a = row lengths TOP->BOTTOM.  Read strip from bottom-left; descent set is
    the set of partial sums of the row lengths read BOTTOM->TOP."""
    b = list(reversed(a))                          # bottom -> top
    D, s = [], 0
    for x in b[:-1]:
        s += x
        D.append(s)
    return sum(D), tuple(D)

def add_strips(mu, e, L):
    """all lambda with lambda/mu a connected e-border-strip, via bead moves."""
    M = maya(mu, L)
    out = []
    for p in sorted(M):
        q = p + e
        if q in M or q >= L + 40:
            continue
        M2 = (M - {p}) | {q}
        lam = from_maya(M2, L)
        ht = len([x for x in M if p < x < q])
        out.append((lam, ht))
    return out

if __name__ == "__main__":
    print(f"{'e':>2} {'mu':>14} {'lambda':>16} {'a(top->bot)':>14} {'r':>2} {'ht':>3} "
          f"{'maj':>4} {'c':>4} {'(r-1)e+maj':>10}  ok")
    bad = 0; total = 0
    rows_for_levelset = []
    for e in (3, 4, 5):
        for mu in parts_upto(9):
            L = len(mu) + e + 3
            for lam, ht in add_strips(mu, e, L):
                r, a = strip_data(mu, lam)
                mj, D = maj_of_strip(a)
                c = n_stat(lam) - n_stat(mu)
                pred = (r-1)*e + mj
                total += 1
                ok = (c == pred)
                if not ok:
                    bad += 1
                    print(f"{e:>2} {str(mu):>14} {str(lam):>16} {str(a):>14} {r:>2} {ht:>3} "
                          f"{mj:>4} {c:>4} {pred:>10}  MISMATCH")
                rows_for_levelset.append((e, r, a, ht, mj, c))
    print(f"\nchecked {total} (mu, strip) pairs for e=3,4,5, |mu|<=9")
    print(f"mismatches against  c = (r-1)e + maj(R):  {bad}")

    # --- is c constant in mu?  level sets ---
    from collections import defaultdict
    byshape = defaultdict(set)
    for e, r, a, ht, mj, c in rows_for_levelset:
        byshape[(e, a)].add(c)
    nonconst = [(k, sorted(v)) for k, v in sorted(byshape.items()) if len(v) > 1]
    print(f"\nribbon shapes whose c is NOT constant over mu: "
          f"{len(nonconst)} / {len(byshape)}")
    for k, v in nonconst[:6]:
        print(f"   e={k[0]} shape {k[1]}: c takes values {v}")

    # --- the level set that IS invariant ---
    byshape2 = defaultdict(set)
    for e, r, a, ht, mj, c in rows_for_levelset:
        byshape2[(e, a)].add(c - (r-1)*e)
    bad2 = [k for k, v in byshape2.items() if len(v) > 1]
    print(f"\nribbon shapes whose  c-(r-1)e  is NOT constant over mu: {len(bad2)}")
    print("\nLEVEL SET  c - (r-1)e  as a function of shape (e=4, all 2^3 shapes):")
    seen = {}
    for e, r, a, ht, mj, c in rows_for_levelset:
        if e == 4:
            seen[a] = (c-(r-1)*e, mj, ht)
    for a in sorted(seen, key=lambda x: (len(x), x)):
        v, mj, ht = seen[a]
        print(f"   shape {str(a):>12}  ht={ht}  c-(r-1)e = {v:>2}   maj = {mj:>2}")
