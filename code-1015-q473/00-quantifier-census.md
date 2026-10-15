# Q473 Task 0 — the quantifier census (written BEFORE any computation)

"The functions lim_{c→∞} G_κ(x_1,…,x_L) are linearly independent" has four unstated
quantifiers. Measured off the source, not guessed.

## (1) Over what field? — DECIDED BY THE DEDUCTION, AND IT KILLS THE BRIEF'S §2.

The relation whose coefficients must be stripped is `eq: G nu integral forumla thm-cinf2`
(line 3857):

    G_ν^∞(x_1..x_L) = Σ_κ A_{ν,κ} · G_κ^∞(x_1..x_L),

    A_{ν,κ} = ∮_{C_1}…∮_{C_n} [∏_{i<j} …] ∏_{i=1}^n [ (w_i−a)(1−qw_i²) / (w_i(1−aw_i)(1−w_i²))
                 · y_{ν_i}/(1−q w_i y_{ν_i}) · ∏_{j=1}^{ν_i−1} (1−w_i y_j)/(1−q w_i y_j) ]
               · lim_{c→∞} F_κ(w_1..w_n).

**A_{ν,κ} depends on the y-alphabet.** (Explicitly through y_{ν_i} and y_1..y_{ν_i−1}, and
again through F_κ.) Therefore the independence required is over a field CONTAINING the y's:

    K := ℚ(q, a, y_1, y_2, y_3, …),   variables = x_1,…,x_L.

⚠ **This invalidates the pre-registered strategy of §2.** Observations A and B of the brief
both separate the G_κ by taking residues/derivatives in the y-alphabet. But y-functionals are
not K-linear — the coefficients c_κ are themselves y-dependent, so Res_{y_m} does not pass
through them. The brief's "filtration by max(κ) is free and I believe it" is TRUE as a
statement about y-support and USELESS as a separating device. The brief's own second route
(§2 last paragraph, "the x-degree/pole structure") is the only admissible one.

Recorded as a brief error, found at census time, before code. This is the
`[[a-sector-is-identified-by-the-parameters-it-carries]]` move: count which symbols the
COEFFICIENT carries before choosing the functional.

## (2) Which index set?
κ ranges over **all finite subsets of ℕ = ℤ_{≥1}** — stated verbatim at line 443: "the sum on
the left is over all finite subsets κ ⊂ ℕ". In `conj:F-orth` the cardinality is further
restricted to |κ| ≤ n = |ν| (because F_κ(w_1..w_n) vanishes for |κ| > n). Writing a subset
decreasingly, κ = (κ_1 > κ_2 > … > κ_d ≥ 1), d = |κ|; κ = ∅ allowed (d = 0).
Open sub-question: is any G_κ^∞ ≡ 0? (outcome O4) — to be measured.

## (3) At which L?
`conj:F-orth` contains **no L**: its LHS is an n-fold contour integral in w_1..w_n only.
L enters solely through the relation in (1), which holds **for every L ≥ n** (hypothesis of
`G nu thm` is L ≥ n). Hence:

    independence at ONE value of L ≥ n suffices to license the coefficient-stripping.

So the theorem to prove is "∃ L" — and I will in fact prove it for **every** L ≥ max|κ|,
which is strictly stronger and costs nothing.

## (4) Finite or infinite family?
The family is **infinite**, and the relation to be stripped is an **infinite** sum (convergent
under `Cauchy identity condition-2`). So there are two distinct statements:

  (LI-fin)  every finite subfamily of {G_κ^∞} is K-linearly independent;
  (LI-top)  every CONVERGENT relation Σ_κ c_κ G_κ^∞ = 0, c_κ ∈ K, has all c_κ = 0.

(LI-top) ⟹ (LI-fin), strictly. The authors' sentence says "linearly independent", which
literally means (LI-fin); their deduction uses (LI-top). **Do not report (LI-fin) as closing
the conjecture without saying which one the deduction needs.** Flagged now so the session
cannot retrofit.
