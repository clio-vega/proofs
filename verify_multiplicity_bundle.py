"""Verify Multiplicity Bundle Theorem - pure Python/NumPy."""
import numpy as np
from itertools import product as iterproduct
from collections import Counter

def make_sigma(n, k, i, j):
    """Build sigma_{[i,j]} on (C^n)^{otimes k}: reversal of positions i..j (1-indexed)."""
    dim = n**k
    sigma = np.zeros((dim, dim), dtype=int)
    tuples = list(iterproduct(range(n), repeat=k))
    idx = {t: m for m, t in enumerate(tuples)}
    for t in tuples:
        t_list = list(t)
        t_list[i-1:j] = reversed(t_list[i-1:j])
        sigma[idx[tuple(t_list)], idx[t]] = 1
    return sigma

def crystal_fi(b, i, n, k):
    """Apply f_i using signature rule for sl_n on a k-tuple."""
    signs = []
    positions = []
    for pos in range(k):
        if b[pos] == i:
            signs.append('+')
            positions.append(pos)
        elif b[pos] == i+1:
            signs.append('-')
            positions.append(pos)
    if not signs:
        return None
    # Parenthesis matching: each - cancels with nearest unmatched + to its LEFT
    matched = [False] * len(signs)
    for idx in range(len(signs)):
        if signs[idx] == '-':
            for q in range(idx-1, -1, -1):
                if signs[q] == '+' and not matched[q]:
                    matched[q] = True
                    matched[idx] = True
                    break
    # f_i acts on LEFTMOST unmatched +
    target = -1
    for idx in range(len(signs)):
        if signs[idx] == '+' and not matched[idx]:
            target = positions[idx]
            break
    if target == -1:
        return None
    result = list(b)
    result[target] = i + 1
    return tuple(result)

def crystal_ei(b, i, n, k):
    """Apply e_i using signature rule for sl_n on a k-tuple."""
    signs = []
    positions = []
    for pos in range(k):
        if b[pos] == i:
            signs.append('+')
            positions.append(pos)
        elif b[pos] == i+1:
            signs.append('-')
            positions.append(pos)
    if not signs:
        return None
    # Same parenthesis matching
    matched = [False] * len(signs)
    for idx in range(len(signs)):
        if signs[idx] == '-':
            for q in range(idx-1, -1, -1):
                if signs[q] == '+' and not matched[q]:
                    matched[q] = True
                    matched[idx] = True
                    break
    # e_i acts on RIGHTMOST unmatched -
    target = -1
    for idx in range(len(signs)-1, -1, -1):
        if signs[idx] == '-' and not matched[idx]:
            target = positions[idx]
            break
    if target == -1:
        return None
    result = list(b)
    result[target] = i
    return tuple(result)

def crystal_decompose(n, k):
    """Decompose B(omega_1)^{otimes k} for sl_n into connected components."""
    elements = list(iterproduct(range(1, n+1), repeat=k))
    hw_elements = []
    for b in elements:
        is_hw = True
        for i in range(1, n):
            if crystal_ei(b, i, n, k) is not None:
                is_hw = False
                break
        if is_hw:
            hw_elements.append(b)

    components = {}
    element_to_comp = {}
    for cidx, hw in enumerate(hw_elements):
        comp = []
        queue = [hw]
        seen = {hw}
        while queue:
            curr = queue.pop(0)
            comp.append(curr)
            element_to_comp[curr] = cidx
            for i in range(1, n):
                nxt = crystal_fi(curr, i, n, k)
                if nxt is not None and nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        components[cidx] = {'hw': hw, 'size': len(comp), 'elements': comp}
    return components, element_to_comp, hw_elements

def count_syt(partition):
    """Count standard Young tableaux via hook length formula."""
    from math import factorial
    n = sum(partition)
    if n == 0:
        return 1
    hooks = []
    for i, row_len in enumerate(partition):
        for j in range(row_len):
            arm = row_len - j - 1
            leg = sum(1 for r in range(i+1, len(partition)) if partition[r] > j)
            hooks.append(arm + leg + 1)
    result = factorial(n)
    for h in hooks:
        result //= h
    return result

def gl_n_dim(partition, n):
    """Dimension of GL_n irrep via hook-content formula."""
    numer = 1
    denom = 1
    for i, row_len in enumerate(partition):
        for j in range(row_len):
            numer *= (n + j - i)
            arm = row_len - j - 1
            leg = sum(1 for r in range(i+1, len(partition)) if partition[r] > j)
            denom *= (arm + leg + 1)
    return numer // denom

def gen_partitions(n, max_length=None):
    """Generate all partitions of n with at most max_length parts."""
    if max_length is None:
        max_length = n
    if n == 0:
        yield ()
        return
    def _gen(remaining, max_part, max_len):
        if remaining == 0:
            yield ()
            return
        if max_len == 0:
            return
        for part in range(min(remaining, max_part), 0, -1):
            for rest in _gen(remaining - part, part, max_len - 1):
                yield (part,) + rest
    yield from _gen(n, n, max_length)

# =====================================================================
print("=" * 70)
print("TEST 1: Crystal decomposition vs Schur-Weyl prediction")
print("=" * 70)

all_pass = True
for n, k in [(2, 3), (2, 4), (3, 3), (4, 3)]:
    print(f"\n--- sl_{n}, k={k} ---")
    components, _, hw_elts = crystal_decompose(n, k)
    crystal_mults = Counter(comp['size'] for comp in components.values())

    sw_mults = Counter()
    total = 0
    for lam in gen_partitions(k, max_length=n):
        d_V = gl_n_dim(lam, n)
        f_lam = count_syt(lam)
        sw_mults[d_V] += f_lam
        total += d_V * f_lam
        print(f"  lambda={lam}: dim(V)={d_V}, f^lambda={f_lam}")

    match = crystal_mults == sw_mults
    print(f"  Crystal mults by size: {dict(sorted(crystal_mults.items()))}")
    print(f"  Schur-Weyl prediction: {dict(sorted(sw_mults.items()))}")
    print(f"  Total: {total} = {n}^{k} = {n**k}: {total == n**k}")
    print(f"  MATCH: {match}")
    if not match:
        all_pass = False

# =====================================================================
print("\n" + "=" * 70)
print("TEST 2: Eigenvalues of sigma_{[1,k]}")
print("=" * 70)

for n, k in [(2, 3), (2, 4), (3, 3)]:
    print(f"\n--- sl_{n}, k={k} ---")
    sigma = make_sigma(n, k, 1, k)
    evals = np.linalg.eigvals(sigma.astype(float))
    eigen_counts = Counter(np.round(evals.real, 8))
    print(f"  Eigenvalues: {dict(sorted(eigen_counts.items()))}")
    print(f"  sigma^2 = I: {np.allclose(sigma @ sigma, np.eye(n**k))}")

    # Commutation with GL_n
    np.random.seed(42)
    g = np.random.randn(n, n)
    g_tensor = g
    for _ in range(k-1):
        g_tensor = np.kron(g_tensor, g)
    print(f"  Commutes with GL_{n}: {np.allclose(sigma @ g_tensor, g_tensor @ sigma, atol=1e-8)}")

# =====================================================================
print("\n" + "=" * 70)
print("TEST 3: Naive reversal vs crystal components")
print("=" * 70)

for n, k in [(2, 3), (2, 4), (3, 3)]:
    print(f"\n--- sl_{n}, k={k} ---")
    components, elem_to_comp, _ = crystal_decompose(n, k)

    mismatches = 0
    examples = []
    for b in iterproduct(range(1, n+1), repeat=k):
        b_rev = tuple(reversed(b))
        if elem_to_comp[b] != elem_to_comp[b_rev]:
            mismatches += 1
            if len(examples) < 3:
                examples.append((b, b_rev, elem_to_comp[b], elem_to_comp[b_rev]))

    if mismatches == 0:
        print(f"  Naive reversal preserves components: YES (all {n**k} elements)")
    else:
        print(f"  Naive reversal preserves components: NO ({mismatches} mismatches)")
        print(f"  (EXPECTED: crystal commutor != naive reversal)")
        for b, br, c1, c2 in examples:
            print(f"    {b} (comp {c1}, size {components[c1]['size']}) -> {br} (comp {c2}, size {components[c2]['size']})")

# =====================================================================
print("\n" + "=" * 70)
print("TEST 4: B(omega_1)^2 multiplicity-free check")
print("=" * 70)

for n in [2, 3, 4, 5]:
    print(f"\n--- sl_{n} ---")
    components, _, _ = crystal_decompose(n, 2)
    sizes = [comp['size'] for comp in components.values()]
    print(f"  Components: {len(components)}, sizes: {sorted(sizes)}")
    print(f"  Multiplicity-free: {len(sizes) == len(set(sizes))}")

print("\n" + "=" * 70)
print(f"ALL TESTS COMPLETE. Overall pass: {all_pass}")
print("=" * 70)
