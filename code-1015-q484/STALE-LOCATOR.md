# STALE LOCATOR in this directory's `q484_engine.py` (md5 370e37b9) -- DO NOT PATCH

Decided in writing, PROVE 2026-10-18 c2.  This directory is a **frozen session
record**: the engine here is the code that actually ran on the date in the
directory name, and the numbers in that date's paper were produced by *this*
file.  Rewriting its docstring would make the archive misreport what ran, which
is a worse fault than a stale locator.  So the correction is recorded alongside,
not inside.

Correction: `motegi_ge`'s docstring in this copy says Motegi's order is defined at
**L640-650** of `2610.15344`.  It is not.  The order is defined at **L582-589**;
L640-650 lies inside the proof of the triangular vanishing and contains no
definition (L670-675 is where that proof concludes).  Checked first-hand against
the e-print on 2026-10-18.

The live copy is `projects/reviews/scratch-20261017c2/q484_engine.py`
(md5 `f53c0258`), which carries the corrected locator *and* the corrected type-B
normaliser `EB`.  **Use that one for new work.**  Copy census over
`/home/clio` + `/tmp` on 2026-10-18 c2: 8 pre-existing copies at 3 distinct
md5s; the wider universe adds nothing to `find projects`, so the enumeration is
now closed.
