"""Littlewood-Richardson coefficients via LR skew tableaux (lattice/ballot words).

c^nu_{lambda,mu} = # of SSYT of shape nu/lambda, content mu, whose reverse
reading word is a lattice word.

Positive controls (untuned, from standard references):
  c^{(3,2,1)}_{(2,1),(2,1)} = 2
  c^{(4,2)}_{(2,1),(2,1)}   = 1
  c^{(2,1)}_{(1),(1,1)}     = 1
  s_{(2,1)}^2 expansion multiplicities (Macdonald I.9 Ex).
"""
from functools import lru_cache
from itertools import product


def _pad(lam, n):
    return tuple(lam) + (0,) * (n - len(lam))


def _trim(lam):
    lam = list(lam)
    while lam and lam[-1] == 0:
        lam.pop()
    return tuple(lam)


def lr_coeff(nu, lam, mu):
    """c^nu_{lam,mu} by explicit enumeration of LR skew tableaux of shape nu/lam."""
    nu, lam, mu = _trim(nu), _trim(lam), _trim(mu)
    if sum(nu) != sum(lam) + sum(mu):
        return 0
    n = len(nu)
    if len(lam) > n:
        return 0
    lam = _pad(lam, n)
    # lam must fit inside nu
    for i in range(n):
        if lam[i] > nu[i]:
            return 0
    if not all(nu[i] >= nu[i + 1] for i in range(n - 1)):
        raise ValueError("nu not a partition")
    if not all(lam[i] >= lam[i + 1] for i in range(n - 1)):
        return 0
    m = len(mu)
    if m == 0:
        return 1 if tuple(nu) == tuple(_trim(lam)) else 0

    # cells of nu/lam, row by row; we fill row by row, left to right,
    # but the ballot condition is on the REVERSE reading word: read rows
    # top-to-bottom, each row RIGHT-to-LEFT.
    rows = [(i, list(range(lam[i], nu[i]))) for i in range(n)]

    count = 0
    # entries[i][j] = value in row i, column j
    grid = [dict() for _ in range(n)]

    def rec(i, j_idx, content, word_counts):
        """Fill row i, position j_idx (index into rows[i][1])."""
        nonlocal count
        if i == n:
            if all(content[k] == mu[k] for k in range(m)):
                count += 1
            return
        cols = rows[i][1]
        if j_idx == len(cols):
            rec(i + 1, 0, content, word_counts)
            return
        j = cols[j_idx]
        for v in range(1, m + 1):
            if content[v - 1] == mu[v - 1]:
                continue
            # weakly increasing along the row
            if j_idx > 0 and grid[i][cols[j_idx - 1]] > v:
                continue
            # strictly increasing down the column
            if i > 0 and j in grid[i - 1] and grid[i - 1][j] >= v:
                continue
            # the entry above (in nu/lam) must exist if row i-1 reaches column j
            # (guaranteed by skew shape); no constraint otherwise.
            grid[i][j] = v
            content[v - 1] += 1
            # ballot condition on the reverse reading word: we generate the
            # reverse reading word in order if we fill each row right-to-left.
            # Here we fill left-to-right, so we check ballot at row end instead.
            rec_ok = True
            if rec_ok:
                rec(i, j_idx + 1, content, word_counts)
            content[v - 1] -= 1
            del grid[i][j]

    # Simpler & safer: enumerate all semistandard fillings, then test ballot.
    cells = [(i, j) for i in range(n) for j in range(lam[i], nu[i])]
    ncells = len(cells)
    if ncells != sum(mu):
        return 0

    result = 0
    filling = {}

    def ballot_ok_partial():
        # reverse reading word: rows top->bottom, each right->left
        cnt = [0] * (m + 1)
        for i in range(n):
            for j in range(nu[i] - 1, lam[i] - 1, -1):
                v = filling.get((i, j))
                if v is None:
                    return True
                cnt[v] += 1
                if v > 1 and cnt[v] > cnt[v - 1]:
                    return False
        return True

    def rec2(idx):
        nonlocal result
        if idx == ncells:
            content = [0] * m
            for v in filling.values():
                content[v - 1] += 1
            if all(content[k] == mu[k] for k in range(m)) and ballot_ok_partial():
                result += 1
            return
        i, j = cells[idx]
        for v in range(1, m + 1):
            if j > lam[i] and filling[(i, j - 1)] > v:
                continue
            if i > 0 and (i - 1, j) in filling and filling[(i - 1, j)] >= v:
                continue
            filling[(i, j)] = v
            if ballot_ok_partial():
                rec2(idx + 1)
            del filling[(i, j)]

    rec2(0)
    return result


if __name__ == "__main__":
    tests = [
        (((3, 2, 1)), (2, 1), (2, 1), 2),
        (((4, 2)), (2, 1), (2, 1), 1),
        (((2, 1)), (1,), (1, 1), 1),
        (((3, 3)), (2, 1), (2, 1), 1),
        (((2, 2, 2)), (2, 1), (2, 1), 1),
        (((4, 1, 1)), (2, 1), (2, 1), 1),
        (((3, 2, 1)), (1, 1), (2, 1, 1), 1),
        (((2, 1)), (2, 1), (), 1),
        (((5, 3)), (3, 1), (2, 2), 1),
    ]
    ok = True
    for nu, lam, mu, want in tests:
        got = lr_coeff(nu, lam, mu)
        flag = "OK " if got == want else "FAIL"
        if got != want:
            ok = False
        print(f"{flag} c^{nu}_{{{lam},{mu}}} = {got}  (expected {want})")
    # s_{21}^2 = s_{42}+s_{411}+s_{33}+2 s_{321}+s_{3111}+s_{222}+s_{2211}
    print("\n s_{21}^2 expansion (Macdonald control):")
    for nu in [(4,2),(4,1,1),(3,3),(3,2,1),(3,1,1,1),(2,2,2),(2,2,1,1)]:
        print("  ", nu, lr_coeff(nu,(2,1),(2,1)))
    print("\nall control tests passed:", ok)
