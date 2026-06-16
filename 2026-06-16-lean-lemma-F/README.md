# Lemma F / Lemma F2 — Lean snapshot (2026-06-16)

Snapshot of `tworow_d4_kernel/TworowD4Kernel/LemmaF.lean` (Lean 4, Mathlib v4.30.0).

Two named declarations, **zero `sorry`, zero warnings, 3 standard axioms**
(`propext`, `Classical.choice`, `Quot.sound`):

- `lemma_F`  — `Q ≥ 4`, single deleted factor, surplus `≥ 1` (the §1.3 factor-in-product engine).
- `lemma_F2` — `Q ≥ 6`, two distinct deleted factors, surplus `≥ 1` (the `c=3` `a`-odd-top keystone;
  `Q ≥ 6` is sharp).

These are the elementary 2-adic primitives behind the three-row `d = 4` boundary lemma
(`proofs/2026-06-15-boundary-lemma-threerow.md` §1.3, `proofs/2026-06-16-c3-boundary-complete.md`).

Honesty note: the LEAN.md stretch wording (`Q ≥ 5`, surplus `≥ 0`) is true but vacuous; the content
form is `Q ≥ 6`, surplus `≥ 1`. Details in `memory/for-robin/2026-06-16-lemma-F-lean.md`.

To build: place under `TworowD4Kernel/`, add `import TworowD4Kernel.LemmaF` to the root, `lake build`.
