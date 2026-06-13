# Three-row c=1 family — d=4 even-|J*| program (2026-06-13)

Verification scripts for `projects/proofs/2026-06-13-threerow-c1-Jstar-even.md`.

- `mn.py`        — Murnaghan–Nakayama (beta-numbers), M_j via characters.
- `dj.py`        — general 3-row closed form M_j = Σ_σ sgn(σ) D_j(...) (JT).
- `c1d.py`       — clean c=1 form M_j = b C(N,b-j+1) - (b-j)C(N,b-j) - (j-1)C(N,b-j-1).
- `prop2b.py`    — factored M_j and Prop-2 val(j) formula (j<=b-1).
- `general.py`   — unified Delta(j)=val(j)-val(0) Kummer formula.
- `finalcheck.py`— compensation lemma + full box + tie conditions (m<=80).
- `checkbdry.py` — boundary families b=1,b=4 hand-formulas; spurious-competitor check.
