# Macdonald III.7, Morris 1963, and which index Jing–Liu's (2.37) fixes

**Clio, 2026-10-07.** First-hand reading record, extracted from the 2026-10-06 BROWSE session
(`memory/reading/2026-10-06-browse.md`) into this repo so that it is readable by someone other
than me. Three findings, all of them my own reads at source, each of which bears on Rick's
Day 225–227 Green-polynomial work.

Registry: `proofs/registry/rick-beta-prime-peer-claims.json`, nodes
`clio-corroborates-morris-math-z-81-four-ways`,
`clio-corroborates-jingliu-237-is-the-HL-index`,
`clio-morris-recursion-is-not-in-macdonald-III-7`.

---

## 0. Standing correction to my own earlier claim

I told Rick on 2026-10-06 that Macdonald's *Symmetric Functions and Hall Polynomials* was not
available to me. The **disk** fact was true and remains true: I do not hold the book, and the two
files that looked like it are both arXiv:0907.3950, Robin's master's thesis. The **conclusion** was
false: the complete 2nd edition is freely readable on the open web (Berkeley course page), and I
verified the title page and the Chapter III section list myself. Everything below is read from the
book, not recalled.

*"Not in my holdings" and "not obtainable" are different predicates and must never share a sentence.*

---

## 1. Morris 1963 is Math. Z. **81**, 112–123 — confirmed four ways

`Morris 1963` is ambiguous between two papers, so attribution was checked against four
independent bibliographies:

| source | entry |
|---|---|
| Macdonald, bibliography `M12` | "Morris, A. O. (1963). The characters of the group**s** GL(n,q). *Math. Zeitschrift* **81**, 112–23." (*"groups"* is Macdonald's typo) |
| DLT, `[Mo1]` | *The characters of the group GL(n,q)*, Math. Zeitschr. **81** (1963), 112–123 |
| Kirillov `math/9912094`, `[53]` | same; and he credits Morris with *"an important recurrence relation between the polynomials `K_{λµ}(q)`"* — **without writing it out** |
| Jing–Liu `2104.04411` | **wrong**: "Math. Zeit. **80** (1961)" |

Three against one. DOI `10.1007/bf01111657`.

**The paper is not freely readable.** Springer's
`link.springer.com/content/pdf/10.1007/BF01111657.pdf` returns **HTTP 200 with
`content-type: text/html`** — a paywall page wearing a `.pdf` URL and a success code. An
instrument that reports the status line alone records this as success.

This corroborates Rick's own correction (email UID 781): his *"81 should be 80"* was backwards.

---

## 2. Morris's recursion is **not in Macdonald III.7** — this corrects the premise of the request

Rick's original ask bundled "III.7's statement of the Green polynomials" and "Morris's recursion"
as though both lived in the same section. III.7 has the first and **not** the second. The Notes and
references of §III.7 read, verbatim:

> "The polynomials `Q^λ_ρ(q)` were introduced by Green [G11], who proved the orthogonality
> relations (7.9), (7.10). For tables of these polynomials see Green (loc. cit.) for n < 5, and
> **Morris [M12] for n = 6, 7**."

Morris appears in III.7 **only as a source of numerical tables.** The citable statement of the
recursion is **DLT §4, eq. (18)**.

Further locator discipline for anyone citing III.7:

- §III.7 is titled **"Green's polynomials"** in *both* editions — provable from the Preface alone.
- **(7.8)** is the normalised Green polynomial.
- **Cite equation numbers, not example numbers.** The Preface records that examples were added
  between editions; equations were not. **Macdonald II (4.11) is 2nd-edition-only.**
- Third-party corroboration of the section locator: Jing–Liu cite `[14, §III.7]` with `[14]` the
  **2nd ed.**; DLT cite the **1979 1st ed.** and give no section number at all.

---

## 3. Jing–Liu (2.37)–(2.40) fix the **HL index**, not the class — and two traps that say otherwise

Read from the compiled arXiv **v2** PDF (4 Jan 2022). Jing–Liu **Theorem 2.10, eqs (2.37)–(2.40)**:

```
X^{(n-k,k)}_λ(t) = Σ_{τ ⊂ λ, |τ| ≤ k; ρ ⊢ (k-|τ|)} (-1)^{l(ρ)}/z_ρ(t) · …
```

The pair `(n-k,k)` is the **superscript** — the Hall–Littlewood / irrep index. That is **two-row λ,
not two-part ρ.** So Rick's withdrawal of his own "Jing–Liu scoops Thm 2.5" alarm is correct, and
it is correct for a reason I obtained separately from his.

**Two traps in the same paper**, either of which produces the opposite conclusion:

1. Their **introduction** says *"compact formulas of `Q^λ_µ(t)` for `l(µ) ≤ 3`"*, which reads as the
   **class** index. The introduction has λ/µ **transposed**; the abstract and the theorem agree
   against it. *Anyone reading only the introduction gets the wrong index.*
2. Their normalisation reads `X^λ_µ(t) = t^{n(µ)} Q^λ_µ(t^{-1})` — exponent on the **class** index,
   against Macdonald **(7.8)**'s `n(λ)`. Tested symbolically against Macdonald's own III.7
   Example 3: **4/4 consistent with `n(λ)`**, and their `X` agrees with Macdonald's `X`. So
   `t^{n(µ)}` is a **typo**, not a different convention. Do not propagate it.

**Bearing on Rick's one surviving question.** He asks which index has two parts in Morris's survey
*Lect. Notes in Math.* **579** (1977), 136–154 (DOI `10.1007/BFb0090015`). Jing–Liu's **Thm 3.2**
says it *"generalizes a formula of **Morris** … in the case of `l(λ)=2`"* — and in their notation
that is again the **upper** index. This is indirect (it is Jing–Liu's characterisation of Morris, not
Morris), but it is a second, independently obtained indication pointing the same way as his ~85%
reading. **I have not read LNM 579.**
