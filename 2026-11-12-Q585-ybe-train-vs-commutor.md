# Q585 at the Schur point — summary (full note: 2026-11-12-Q585-ybe-train-vs-commutor.tex/.pdf)

**Verdict: ill-posed at the instrument.** The train built from the KT puzzle model's own crossing is not deterministic, so the falsifier ("c≥2 case with T≠C, T deterministic") cannot fire. The prediction is neither confirmed nor refuted.

- **Prop (proved):** the Schur-point SE×SW crossing D0 = ∇∘Δ factors through one 3-label edge, so its rank is exactly 3 of 9. It is not a braiding.
- **Lemma (proved):** D0 kills the input (SE,SW)=(0,1), which is the equivariant-rhombus configuration, whose weight vanishes at y=0.
- **Claim (computed, N=2..5):** the bundle crossing B_N is never single-output, and it realises the swap only on the constant words.
- **Commutor C (computed):** two instruments agree 1308/1308 (n=4, |λ|,|μ|≤5, 67 cells with c≥2). The lesion differs on 1298. The lex control differs on 67/67 cells; reverse-lex agrees on 67/67, so it would have been blind.
- **Gap / next:** a train with a *foreign* defect R is not ruled out. Defining one needs a source that states the train argument. The candidate is the HK octahedron commutor on hives (`math/0408114`, recalled, unverified).
- Code: `code-1112-Q585/` @ `4641a0145b29`. Registry: `registry/q585-ybe-train-vs-commutor.json`.
