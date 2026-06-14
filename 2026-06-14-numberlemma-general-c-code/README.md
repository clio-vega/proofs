# NL_c verification scripts (2026-06-14)

Number Lemma NL_c:  v2 C(F+2c-1,j) + v2(j^(2c)) >= v2(F) + beta(c),  beta(c)=(c-1)+v2((c-1)!).

- `nlc_verify.py`   — beta closed form; NL_c + sharpness (min slack 0), c=1..5.
- `nlc_steps.py`    — subset identity; anchor identity = sum_{i=0}^{2c-1} v2(F+i); min over F = beta(c).
- `s1_test.py`      — Gap1 T1/T2 split both >= v2(b+3)+1 (a even); j<6 region.
- `s1_termwise.py`  — term-by-term subset-identity bound for the heavy part C(b+3,j)H(j).

All: 0 failures. Run `python3 <name>.py`.
