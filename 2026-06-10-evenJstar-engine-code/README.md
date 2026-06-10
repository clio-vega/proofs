# Reproducing scripts — 2026-06-10 π-adic engine for the d=4 fiber

Pure-Python (sympy for one symbolic check; Murnaghan–Nakayama in `symfunc.py`/`job1_tie_census.py`).
Run any `job_*.py` from this directory. See `../2026-06-10-evenJstar-box-steplaw.md` §6 for the
claim→result table. Key scripts:

- `job_mod2_engine.py` — the engine `Φ(z)≡χ^λ(2^m)(1+z)^m mod 2` and the box reduction (Test 1/2).
- `job_verify_engine.py` — exact `(1+z)`-lift `Φ=Σ C(m,r)2^r R_r(1+z)^{m-r}`, `R_0=χ^λ(2^m)`.
- `job_sharp_mod2.py` — sharp `J*` bottom-mod-2 `= y^{l0}(1+y)^{g}` (1624/1624 ties).
- `job_hierarchy.py` — exact `v_π(G)`; coarse vs sharp polygon minima.
- `job_bpicture.py`, `job_submask.py`, `job_box_detail.py`, `job_v2D_explore.py` — exploration.
