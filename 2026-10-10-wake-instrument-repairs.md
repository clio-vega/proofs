# WAKE 2026-10-10 — two instrument repairs, and what they were hiding

Written because `projects/code/` and `projects/memory/` are **not git repositories**
(`git rev-parse` exits 128 in both), so neither repair is readable by Robin or Rick
where it was made. This note lives in `clio-vega/proofs` so the findings survive.

---

## 1. `citation_check.py` reported 176 false problems from a stale constant

`code/citation_check.py:21` read

```python
LEVELS = ["agent-summary", "abstract", "deep-read", "verified-quote"]
```

while **176 entries** in `memory/reading/sources.json` carried `extraction: "title-only"`.
Every run therefore emitted 176 `extraction 'title-only' invalid` problems.

`title-only` is **deliberate**, and `sources.json`'s own `schema_notes` says so
(`"2026-10-04-c5"`): wrong-ID near-misses and planted control IDs *"are now registered at
title-only with their ROLE stated in notes"*. For a negative record it is the **only honest
value** — there is no abstract I read. **The constant was stale, not the corpus.**

Fixed by adding `"title-only"` to `LEVELS`, with the reasoning in a comment so a future
session does not "fix" the 176 entries instead.

**Readings.** One file (`topics/methodological-principles.md`):

| | problems | extraction complaints |
|---|---|---|
| before | 436 | 177 |
| after  | **255** | **0** |

`436 = 255 + 176 + 5` — the 176 stale false positives plus the 5 orphans registered below.
Two-file arm (adding `dream-journal/2026-10-09.md`): **442 → 256**, a drop of 186 = 176 + **10**,
the 5 orphans counted **twice** because they appear in both files. That double-count is exactly
what the 10-09 dream observed and could not explain.

**Positive control.** Planting `extraction: "skimmed"` is still flagged, and the message now
prints the widened list — so the enum was **widened, not disabled**. `sources.json` md5-restored
byte-identical after the control arm.

**Consequence worth stating.** Every `citation_check` reading I have taken since `title-only`
came into use has been inflated by 176. The "years-deep citation backlog" I have been reporting
and declining to touch is real but **substantially smaller** than the numbers I quoted for it.
A stale constant in my own validator was manufacturing most of the backlog.

## 2. Five arXiv IDs lived in prose since May 2026 with no index entry

`0709.3766`, `2305.01306`, `2504.12798`, `2507.10061`, `2605.17844` — reported every cycle,
**deferred four times with the diagnosis recorded in place of the action.**

The reason the deferral kept happening is real: an honest `extraction` level needs each context
read, and **the contexts record the conclusion ("deflated"), not the reading depth.** The levels
had to be reconstructed from the May browse logs, not the incident notes. Four are
`agent-summary`; one is `title-only`. **None is `deep-read` — none was read first-hand.**

All are incident records from a **different domain** (Hecke-categorical / W-graph), not sources
in my current LR/cylindric territory: four are the worked instances behind four methodological
principles, and `0709.3766` is a **negative record** for an ID I once hallucinated into my own
agent prompt (it is a Phys. Rev. A entanglement paper; the Zinn-Justin seed is `0809.2392`).

**Titles were fetched from arXiv, not copied from my own prose** — deliberately, because a
paraphrase in a title field renames the object. That caught two errors in my own records:

- `2605.17844` is **27 pages**, not the *"22-page paper"* asserted in two of my files — and the
  page count was the thing I had cited as evidence the paper was real.
- `2504.12798` and `2305.01306` are **Ho *and* Li**; my notes carry both under *"Ho"* alone.

**And one error I made while writing the fix.** My first version of the `0709.3766` entry gave
read path `dream-journal/2026-08-29-dream.md` — **a filename I composed rather than checked.**
The real locations are `for-robin/2026-08-29-dream.md` and `dream-journal/2026-09-17-c2.md`.
`citation_check` caught it within minutes. In the one entry whose entire purpose is to record
that I once invented a plausible arXiv **ID**, I invented a plausible **file path**. It is logged
in that entry's `corrections`.

Index: **808 → 813** sources. Verified by an independent re-read in a fresh process, not by the
writing script's own report.

---

## Named ask to Robin — one line I cannot repair myself

`/home/clio/scripts/boot-prompt.md:77` prescribes

```
--files-dir proofs
```

for `trustcheck validate`. That resolves to `proofs/proofs/…` and reported **69 phantom
problems** this morning on a registry that is in fact clean (`--files-dir .` → `0 problems`,
positive control fires). The registry's `file` fields already carry their own `proofs/`,
`lean/` and `reviews/` prefixes, so the correct value is `.`.

`scripts/lean-prompt.md:67–68` and `scripts/peer-review-prompt.md:77–79` **already carry the
correction** — the latter even tabulates `--files-dir proofs` as WRONG. The only prompt still
carrying the defect is the one I cannot edit: `boot-prompt.md` is a **read-only bind mount**
(`open('r+')` → `OSError errno 30`). WAKE is the one phase whose own prompt it cannot repair,
which is why this correction has been rediscovered rather than applied.

**Ask:** change `boot-prompt.md:77` from `--files-dir proofs` to `--files-dir .`.
