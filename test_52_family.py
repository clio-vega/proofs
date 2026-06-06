"""Pin down the (k, 2, 1^q) pattern."""
import sys
sys.path.insert(0, '/home/clio/projects/scratch/2026-05-13-multi-corner')
from compute_dim_B import compute_dim_B, fmt_lam

shapes = [
    ((5, 2), 7),
    ((6, 2), 8),
    ((6, 2, 1, 1), 10),
    # Also see if there's a pattern for (k, k, 1^q)
    ((5, 5), 10),  # too big maybe
]
for lam, n in shapes:
    if sum(lam) != n: continue
    if n > 12: continue
    try:
        dim_B, d_lam = compute_dim_B(lam, n, verbose=False)
        print(f"{fmt_lam(lam):14} n={n:2} dimV={d_lam:5} dim_B={dim_B}", flush=True)
    except Exception as e:
        print(f"{fmt_lam(lam)} ERROR {e}", flush=True)
