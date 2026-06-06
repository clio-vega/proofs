"""Extra tests to refine the paper: q=0 shapes and (4,4,1^q) at higher q."""
import sys
sys.path.insert(0, '/home/clio/projects/scratch/2026-05-13-multi-corner')
from compute_dim_B import compute_dim_B, fmt_lam

def f_lambda(lam):
    from math import factorial
    if not lam: return 1
    n = sum(lam)
    prod = 1
    for r, lr in enumerate(lam):
        for c in range(lr):
            arm = lr - c - 1
            leg = 0
            for rp in range(r+1, len(lam)):
                if lam[rp] > c:
                    leg += 1
                else:
                    break
            prod *= arm + leg + 1
    return factorial(n) // prod

shapes = [
    ((3, 3), 6),
    ((4, 2), 6),
    # bad-pair test for (4,2,2,1)
    ((3, 2, 1, 1), 7),  # = mu_{1,2}
    ((4, 2, 1, 1), 8),  # for cross-check, already stable
    # Test (5, 2, 1, 1) - is this stable or not?
    ((5, 2, 1, 1), 9),
    # And (5, 2, 1, 1, 1)
    ((5, 2, 1, 1, 1), 10),
    # what about (4, 4, 1^q) at larger q to find stability?
    # (4,4,1^4) is too big (n=12, dim 1485). Try (4,4,1^4) timing
]
print(f"{'lambda':18} {'n':>3} {'r_lpc':>5} {'q':>3} {'dimV':>6} {'f(>=2)':>7} {'dim_B':>6}", flush=True)
for lam, n in shapes:
    if sum(lam) != n: continue
    lam_ge2 = tuple(b-1 for b in lam if b >= 2)
    f_pred = f_lambda(lam_ge2)
    r_lpc = 0
    ell = len(lam)
    for i in range(ell - 1):
        if lam[i] > lam[i+1]:
            r_lpc += 1
    q = sum(1 for x in lam if x == 1)
    try:
        dim_B, d_lam = compute_dim_B(lam, n, verbose=False)
        print(f"{fmt_lam(lam):18} {n:>3} {r_lpc:>5} {q:>3} {d_lam:>6} {f_pred:>7} {dim_B:>6}", flush=True)
    except Exception as e:
        print(f"{fmt_lam(lam):18} ERROR {e}", flush=True)
