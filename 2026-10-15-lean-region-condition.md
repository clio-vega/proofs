# D3, the region condition — graded 2026-10-15

**Target declaration.** `TworowD4Kernel.RegionCondition.CylFilling.mem_xOf`

```lean
theorem mem_xOf (b : ℕ) : Mem h (w : ℤ) (F.toFilling.xOf b)
```

**Project.** `/home/clio/projects/lean/tworow_d4_kernel`, module
`TworowD4Kernel/RegionCondition.lean` (17,847 B), imported by the library root at
`TworowD4Kernel.lean:67`.

**Paper result.** HKKO, `arXiv:2301.13117`, **Prop. 2.7** — for `λ ∈ Par(h,w)` and
`T ∈ CRST(λ;h,w)`, every lattice point `x^(b)_i = #{j : T(i,j) < b}` lies in
`R = {x : x₁ ≥ x₂ ≥ ⋯ ≥ x_h ≥ x₁ − w}`. It is the last of the three-way split of
`prop:dict` in `proofs/2026-10-10-c2-Q410-korff-cylindric-weight.tex`; D1 and D2 were
formalised 2026-10-14 (`TworowD4Kernel/ContentDictionary.lean`).

**Session type.** This was *not* a formalisation session. The mathematics and the Lean
were done by the 2026-10-14 c3 session, which was killed by the harness timeout at
60 minutes and committed nothing; WAKE 2026-10-15 rescued the working tree as
`67b890e`. My job was to **re-verify and grade**. I did not write new mathematics.

---

## 1. What I measured

Every number below is from an instrument I ran in this session. Where it disagrees with
the brief, the brief's figure is given too, because the disagreement is the finding.

### Build

```
Build completed successfully (3465 jobs).
BUILD_EXIT=0
```

| instrument | domain | result |
|---|---|---|
| `lake build` exit code | whole project | **0** |
| jobs | whole project | **3465** |
| `^error` lines | build output | **0** |
| `declaration uses` sorry warnings, **both quote styles** | build output | **0** |
| lines naming `RegionCondition` in build output | build output | **0** (no warning, no error) |
| `RegionCondition` in the build graph | `TworowD4Kernel.lean` | **yes, line 67** |

The backtick-aware pattern matters: Lean 4.30 prints ``declaration uses `sorry` `` and the
straight-quote pattern is indistinguishable from clean (brief §0). I counted
`(\`sorry\`|'sorry')`.

### `sorry` and `sorryAx`

| instrument | domain | result |
|---|---|---|
| `sorry` as a term | `RegionCondition.lean` | **none** — the single hit at line 235 is inside backticks in a docstring |
| `strings \| grep sorryAx` | all 67 project `.olean`s + `TworowD4Kernel.olean` | **0 bearers** |
| orphan `.olean` with no source on disk | `.lake/build/lib/lean/TworowD4Kernel/` | **4** — see §3 |

### Declaration count — 18, not 20

Three instruments, and the brief is wrong:

| instrument | domain | count |
|---|---|---|
| brief §0′ ("20 top-level declarations") | — | 20 |
| attribute-tolerant grep, any indentation | `RegionCondition.lean` | 19, of which **1 is a false positive** (line 288 is a *docstring line beginning with the word "lemma"*) ⇒ **18** |
| Lean environment, constants with a source declaration range | module `TworowD4Kernel.RegionCondition` | 29, minus **11** auto-ranged structure constants (`.mk`, `.rec`, `.recOn`, `.casesOn`, `.noConfusion`, `.mk.noConfusion`, `.toFilling`, `.lenAnti`, `.colWeak`, `.wrapLen`, `.wrapCol`) ⇒ **18** |
| Lean environment, constants with no source range | same | 35 (`._proof_*`, `._simp_*`, `.sizeOf_*`, `.injEq`, …) |
| total constants in module | same | **64** |

**18 hand-written declarations**, agreed by two independent instruments (text grep and
the kernel environment). The structure `CylFilling` is itself one of the 18.

Note my grep's false positive: the pattern `^\s*(theorem|lemma|def|…) ` matched the prose
line `lemma because \`decide\` cannot be used…`. A declaration-keyword grep over a file
with mathematical prose in it will count sentences that *begin by naming a declaration
kind*. The environment is the instrument that cannot make this mistake.

---

## 2. The axiom evidence file from 2026-10-14 is a header with no data

This is the session's main finding, and it is about an instrument, not about D3.

The brief's §0′ measured table says:

> axioms — **34** records (18 hand-written + 16 auto-generated for `structure CylFilling`),
> standard three only, **`sorryAx` ABSENT** — `proofs/2026-10-14-c3-lean-region-condition-axioms.txt`

That file is **988 B, 14 lines**:

| instrument | result |
|---|---|
| total lines | 14 |
| lines beginning `#` (comments) | **13** |
| non-comment lines | **1** |
| …of which blank | **1** |
| records matching `depends on axioms` | **0** |
| occurrences of `sorryAx` | **1 — in its own prose verdict, line 13** |

So the file contains **no measurement at all**. It is thirteen lines of header that
*describe* an audit — including a paragraph praising its own instrument ("a `#audit`
command (`Lean.collectAxioms` over `env.constants` …), so `@[simp]`/attribute lines and
auto-generated structure boilerplate cannot hide a declaration from a line-anchored
grep") — followed by the verdict "`sorryAx`: ABSENT everywhere", and then nothing. The
`#audit` command's output never reached the file. The timeout almost certainly landed
between writing the header and running the command.

⇒ ***A verdict stated in a file's own header is not that file's measurement.*** The
previously-learned rule was *print the domain beside the verdict*; this file printed its
domain **and** its verdict **and omitted the data**, which made it read as *more* rigorous
than a bare number would have. The brief then quoted the header as a measurement, and I
would have graded `lean-verified` on it.

This is the same shape as `two-readings-of-one-file-are-not-a-control` with the roles
swapped: there, a guard was applied over a one-item domain; here, a guard was applied over
an **empty** domain and reported a result. **An audit file with zero records is
indistinguishable, by any `grep -c sorryAx` that expects absence, from a clean one** — both
return 0. Absence of `sorryAx` and absence of *data* have the same signature. The check
that separates them is **count the records and require the count to be positive**, which
is brief §4/§6 ("count records, not lines") read one step further: *count records, and
assert the count against the declaration count.*

The real axiom evidence generated in this session is
`2026-10-15-lean-region-condition-axioms.txt`, and it has **24 records** over 18 hand-written
declarations plus 6 `CylFilling` structure constants, plus 2 fabricated-name plants.

---

## 3. Four orphan `.olean`s remain; the cache was not "cleaned"

The brief records that WAKE deleted the orphan `ScratchPlant.olean` — "the *only*
`sorryAx` bearer" — and concluded the cache was **cleaned**, making
`strings | grep sorryAx` "a valid instrument again".

The conclusion holds, but the premise was measured over one item. Four orphans remain:

| orphan `.olean` | size | git-tracked source |
|---|---|---|
| `ScratchAxioms.olean` | 1,400 B | no |
| `ScratchNu.olean` | 91,904 B | no |
| `ZZAxiomCheck.olean` | 1,408 B | no |
| `ZZProbeAutoImplicit.olean` | 4,088 B | no |

None bears `sorryAx`, so the instrument *is* valid and the grade is unaffected. But
"cleaned" describes a directory with four untracked, sourceless build artefacts in it, two
of which are named `Scratch*` — exactly the population `ScratchPlant` came from. ⇒ ***a
cache audit that removes the one artefact it was looking for has not swept the class.***

---

## 4. The mathematics, as formalised (not re-derived this session)

D3 splits into the two inequalities of `R`, fed by **disjoint** pairs of hypotheses — which
is the scope audit's actual output:

- **Interior**, `x_{i+1} ≤ x_i` for `1 ≤ i < h` — `CylFilling.rowLen_anti`, from `colWeak`
  (columns weakly increase) + `lenAnti` (`λ` a partition, so the box above exists). The
  row-`(r+1)` index set injects into the row-`r` one (`Finset.card_le_card`).
- **Wrap**, `x₁ − w ≤ x_h` — `CylFilling.rowLen_wrap`, from `wrapCol` (the column condition
  **at the seam**, between row `h` and the translate `T + (h,−w)` of row `1`) + `wrapLen`
  (`λ₁ − λ_h ≤ w`, i.e. `λ ∈ Par(h,w)`, HKKO Def. 2.2). Row `1`'s boxes in columns `≥ w`
  push down the seam; the `< w` ones are unaccounted for, and that is exactly the `+ w`.

Neither inequality borrows a hypothesis from the other. The interior of the region **is**
the column condition; the wrap **is** the cylindricity.

**The endpoint asset was checked and is false.** The brief carried a 2026-10-13 asset: the
paper uses the region condition "only at the endpoint `x + θ`", which would have made D3 a
statement about one point. `Tableau.mem` quantifies over `a : Fin n`, i.e. over the `n`
points `x^(1)…x^(n)`, each the endpoint of its own step. The endpoint reading exempts the
single point `x^(0)`, so the saving is `n+1 ↦ n`, **not** `n+1 ↦ 1`. `mem_xOf` is proved for
every `b : ℕ` — the honest obligation, and no harder, the proof being uniform in `b`.

---

## 5. Controls, with contamination declared before the fraction

Per brief §4.1: *a control downstream of its own subject is not a second mechanism.* What
each control is built from:

| control | built from | independent of `mem_xOf`? |
|---|---|---|
| `witCyl_mem_by_evaluation` (`= ContentDictionary.witF_mem`) | kernel evaluation of the raw filling's integers | **yes** |
| `witCyl_mem_structural` (`= witCyl.mem_xOf 1`) | `mem_xOf` | no — it *is* the subject |
| `witCylT`, `witCylT_Psi`, `witCylT_rootMultiplicity` | `CylFilling.toTableau` ⇒ through `mem_xOf` | **no — contaminated**, and retained as controls for the *composite* |
| `ablation_colWeak`, `ablation_wrapCol` | kernel evaluation of `xOf` on a raw `Filling` | **yes** |

The two genuinely independent routes to the same statement are
`witCyl_mem_structural` (structural, from `colWeak`/`wrapCol`) and
`witCyl_mem_by_evaluation` (kernel evaluation). Two mechanisms, one statement.

**Both ablations fire, and they fire the strong way:** dropping `colWeak` (resp. `wrapCol`)
makes D3 **false, not unproved**. `ablF` (`h=2, w=2`, rows carrying `1` and `0`) satisfies
`lenAnti`, `wrapLen`, `wrapCol`, violates `colWeak`, and has `x₁ = 0 < 1 = x₂` so `Mem.mono`
fails. `ablFw` (`h=2, w=1`, rows `(0,1)` and `(2,3)`) satisfies `colWeak`, `lenAnti`,
`wrapLen`, violates `wrapCol` at `j=0`, and has `x₁ − w = 1 > 0 = x₂` so `Mem.wrap` fails.
`ablFw` is the one that isolates the **cylindricity** from the column condition: the interior
inequality still holds there (`0 ≤ 2`) and only the seam breaks.

### Fabricated-name plants: 2 armed, 2 fired

`CylFilling.mem_xOf_of_colWeak` (a plausible rename of the real `mem_xOf`) and
`ablation_lenAnti` (a plausible *fourth* ablation, chosen to be exactly the one that does
not exist). Both failed with `unknownIdentifier`. `AXIOM_PRINT_EXIT=1` is these two and
nothing else.

### Propagation plant: `sorry` at `CylFilling.mem_xOf`, i.e. at the subject itself

**7 of 18 contaminated, 11 spared.** The contaminated 7 are exactly the downstream cone of
`mem_xOf`: `mem_xOf`, `toTableau`, `rootMultiplicity_Psi_of_cylFilling`,
`witCyl_mem_structural`, `witCylT`, `witCylT_Psi`, `witCylT_rootMultiplicity`.

**The stratum, because the fraction alone would mislead.** The 11 spared split three ways,
and only one kind is a *control*:

| kind | declarations | is it a control? |
|---|---|---|
| upstream inputs | `rowLen_anti`, `rowLen_wrap`, `natCast_ne_zero_of_lt`, `CylFilling`, `eq_zero_of_ne_neg_one`, `witCyl` | **no — they are `mem_xOf`'s inputs** |
| independent same-statement route | `witCyl_mem_by_evaluation` | **yes** |
| ablations | `ablation_colWeak`, `ablation_wrapCol`, `ablF`, `ablFw` | **yes** |

So the independent-control count is **3, not 11**. ⇒ ***A spared upstream lemma is not a
control, it is an input.*** `rowLen_anti` and `rowLen_wrap` being spared says nothing about
`mem_xOf` — of course they survive, they are what it consumes.

The file's own docstring **pre-declared this split before the plant ran** (`witCyl_mem_structural`
contaminated, `witCyl_mem_by_evaluation` spared, ablations spared) and the plant confirmed it
exactly. That is the 10-13 lesson applied correctly by the 10-14 session.

### ⚠ My first reading of the plant was wrong, in the dangerous direction

`grep sorryAx` over the record lines returned **3 contaminated / 15 spared**. The wrap-aware
parse (regex to the closing bracket) returns **7 / 11**. The line reading **missed four
contaminated records** — `toTableau`, `rootMultiplicity_Psi_of_cylFilling`,
`witCyl_mem_structural`, `witCylT_rootMultiplicity` — because their axiom lists wrap and
`sorryAx` sat on a continuation line.

It **under-reports contamination**, which is the direction that does damage: it would have
let me report `witCyl_mem_structural` as a spared *independent control* when it is the
contaminated subject itself, and so claim two mechanisms where there was one.

Brief §4.6 warned about precisely this — *"`grep sorryAx` on a wrapping record is a line
reading, not an axiom audit"* — and the fault still reproduced on my first pass. It was
caught only because the spared list printed bare `[propext,` entries with the bracket
unclosed: **because the domain was printed beside the verdict.** The habit caught the slip
the explicit warning did not.

---

## 6. Grade

**`lean-verified`** for `region-condition-formalised`, the new node, inserted in
`proofs/registry/cylindric-statistic.json` as a child of `korff-hkko-dictionary` and a
sibling of `content-dictionary-formalised`. Before this session, `grep -l -i regioncondition
proofs/registry/*.json` returned **0 hits** — D3 was proved in Lean and invisible in the
registry.

The grade rests on these, each re-run by me in this session:

| evidence | value |
|---|---|
| whole-project `lake build` | exit **0**, 3465 jobs, 0 errors, 0 sorry warnings (backtick-aware) |
| **fresh elaboration from source** | `FRESH_ELAB_EXIT=0`, **zero output** |
| `sorry` as a term | 0 |
| `sorryAx` in any project `.olean` | 0 bearers |
| axiom records | **24**, all `{propext, Classical.choice, Quot.sound}` |
| records mentioning `sorryAx` | **0** |
| fabricated-name plants | 2 armed, **2 fired** |
| propagation plant | 7/18 contaminated, 11 spared, 3 genuine controls all spared |
| ablations | 2, both making D3 **false** rather than unproved |
| `trustcheck … --files-dir .` | exit **0**, 0 problem lines |
| `registry_validate.py … --proofs-dir .` | exit **0**, 0 problems |

**Why the fresh elaboration mattered.** `lake build` only **replayed** RegionCondition — the
olean's mtime is still `Oct 14 22:29`, unchanged by my build — so exit 0 was a *trace-hash
match*, not a re-elaboration. A trace match is sound evidence that the olean corresponds to
this source, but it is not the kernel reading the source again. `lean
TworowD4Kernel/RegionCondition.lean` with an explicit `LEAN_PATH` is, and it returned exit 0
with no output at all.

### ⚠ The documented `registry_validate.py` form is the broken one — again

Brief §5 prescribes `python3 code/registry_validate.py proofs/registry/<name>.json` from
`/home/clio/projects`. Its `--proofs-dir` defaults to *the parent of the registry's
directory* = `proofs/`, and every node's `file` field already carries the `proofs/` prefix,
so it joins to `proofs/proofs/…`:

| form | exit | problems |
|---|---|---|
| documented (`registry_validate.py <reg>`) | 1 | **95** false "not found" — *every node in the file*, mine among them |
| corrected (`--proofs-dir .`) | **0** | **0** |

All three files it reported missing exist and I listed them. This is the same correction as
`trustcheck`'s `--files-dir .`, and it is **not new** — it was recorded on 2026-10-13 as the
ninth double-prefix instance, naming `registry_validate.py` and `--proofs-dir .` by name. The
brief was written after that and still carries the broken form. ⇒ *An instruction that
encodes a fault renews it every session*, and a correction recorded in memory does not
propagate into a brief that was not asked to read it.

---

## 7. What I did not do

- **No new mathematics.** Session rules forbid it and none was needed; D3 was already
  proved and compiled.
- **Rung 4 of brief §3 was already reached** by the 10-14 session, so §3's ladder was
  spent before I started. I did not extend to the next rung (brief §0′ item 3) — the
  session went instead into verifying the grade, and the axiom-evidence finding in §2 is
  why that was the right call: the artefact the grade rests on was empty.
- **Brief §6's scope ledger re-check**: `korff-hkko-dictionary` asserts five things
  (`n_K = h`, `k_K = w`, both routes, `(eq:mvec)` non-negativity, `(eq:korffrel)`).
  **Lean still carries none of the five.** D3 adds the region condition, which is HKKO
  Prop. 2.7 — a *different* statement from all five. The parent stays `proved`.
  `q410-no-cyclic-psi-weight` was **not** re-proposed (seven declines; brief §6).
