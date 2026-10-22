# For Robin — 2026-10-21 c2: Theorem A audit, one citation fix (and a validator fault)

Short note, because the mathematics did not change.

## The slot was issued twice

This LEAN session was handed a brief whose target c1 had already finished. Armed 00:15,
target committed and pushed **09:19**, this session started **09:54** with the brief
still live. I did not re-formalise it — I audited it instead, which is the only verb
left that adds information. The brief is now archived so it cannot fire a third time.

## Theorem A holds up

`TworowD4Kernel.CFactorisation.T_c_eq` is sorry-free and I confirmed it independently:
`lake build` green (928 jobs), **41** constants in the namespace, **none** depending on
anything outside `{propext, Classical.choice, Quot.sound}`, and the declaration census
closes as a set in both directions (all 28 source declarations present; the 13 extras
all auto-generated).

The one thing I want to flag as *method*: c1's first sorry-detector had given a false
all-clear on a planted sorry, so I did not trust its replacement on its own say-so
either. I re-planted `theorem planted_control (n : Nat) : n = n := by sorry` and
confirmed the detector fires (`[sorryAx]`) on the plant and stays silent on `T_c_eq`.
An all-clear from an untested detector is worth nothing; that control is the only reason
this one is worth something.

Scope is unchanged, and I checked c1's banner against the source rather than trusting
it: `G` and the `q`-power really are abstract, and **`T^1 = δ` really is an input. This
does not prove G5.**

## The fix that matters to you

The module cited **`arXiv:2409.06217`** for Motegi. That is **not his paper** — the
correct ID is **`arXiv:2610.15344`** (*Type B- and C- trigonometric weight functions*).
It was already pushed, so it was the version on GitHub, i.e. the version you'd read.

Fixed in `9f0a327`:
https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/CFactorisation.lean

Three independent confirmations before I touched it: the paper proof's own citation
(`2026-10-20-c2-G5-reduction.tex` L63-64), a corpus count (**580** occurrences of the
right ID across my record, **0** of the wrong one outside that single line), and my
sources index (`authors: "Motegi, Kohei"`). I also sharpened the locators while I was
there — `prop:allcard-Bhat-inverse` statement **L4896** vs proof head **L4910**.

Why it survived the session that wrote it: **everything around it was right.** L4666 is
the correct locator for `lem:allcard-Bhat-product`. And a well-formed wrong arXiv ID
resolves to *some* real paper, so a reader following it lands somewhere plausible
instead of on an error page. The generalisation I'm keeping: a citation slot is graded
by its neighbours, and its neighbours cannot grade it. The cheap instrument is the
corpus count — the right ID scored 580, the wrong one scored 1.

## One thing I found and deliberately did NOT fix

`registry_validate.py`'s source check is **vacuous** on `motegi-typeBC-multiplication.json`.
It prints `OK: ... is valid` with 21 warnings of one shape.

I first guessed a loader bug and **that guess was wrong** — line 123 unwraps the nested
`sources` key correctly. The real cause is a format mismatch: the registry writes
`sources` as ID-plus-locator prose (`"arXiv:2610.15344 eq:allcard-Bhat-product-C L4668, …"`)
while the index is keyed on the bare `"2610.15344"`, so no lookup ever resolves. Positive
control: rewriting them to bare IDs takes that warning class from **21 → 0**.

The consequence is not cosmetic. `known` at line 293 stays empty, so the README's own
boundary rule — *a `proved` node whose cited sources are all `agent-summary` gets
flagged* — **cannot fire on any node of that file**.

**Severity, measured rather than assumed:** all three cited IDs are at `verified-quote`,
so the rule would not have fired anyway. It is **latent, not live** — nothing is
currently hiding behind it. I left it alone because the repair is a format decision on a
shared instrument (split the field, or teach the validator to parse the leading ID), and
that is a review-slot call, not something to make unilaterally from a Lean slot. Your
mirror registry may have the same shape — worth a look.

Artifact: `proofs/2026-10-21-c2-lean-c-factorisation-audit.md`.
