# Q585 G1 — Exchange proved from two two-move identities (PROVE 2026-11-15)

Full write-up: `proofs/2026-11-15-Q585-G1-mixed-commutation-exchange.tex` (+pdf). Code: `proofs/code-1115-Q585-G1/` @ `e10161532c26`.

## Theorem (conditional; proved)
Assume, for all engine states σ:
- **(MC)** for distinct outer corners x≠y of the A-nest: σ;x→B;y→C = σ;y→C;x→B (same tiling, same landing boxes);
- **(SW)** for an outer corner x=(c,r), r≥1, with y=(c,r−1) an outer corner after x leaves: σ;x→B;y→C = σ;x→C;y→B, and the C-landings agree.

Then (Exchange) holds: the first m−α_ℓ landings of (iii) on M equal the landings of (iii) on M_1 (M_1 = M after (ii) sends the top row).
Consequently, by induction on ℓ(α) and the c1 forcing lemma, R_iii = ι_* R_ii (G1, Purbhoo Prop 5.1 (ii)=(iii), in the engine).

## Proof outline
1. Lemma orders: φ(c,r) = (c,r+1) for c<α_ℓ (identity otherwise) preserves the (iii) sort key (depth, −c). So D = φ(D̂) followed by the bottoms (c,0), c<α_ℓ.
2. Chain invariant: σ'_k = σ_k ; Z_k→B (east to west), where Z_k = the tops of the columns c<α_ℓ.
   - Case (a), a cell in a column ≥ α_ℓ: commute it past Z→B using MC.
   - Case (b), x = z_c: MC past Z^W, then SW at σ_k;Z^E→B, then MC past Z^E.
   - All corner conditions follow from both 𝔞_k and 𝔞_k∖Z_k being ideals, plus the cornerness of the M_1 move y.
3. Exchange is the landing clause of the chain. G1 follows by induction: values < ℓ come from the induction hypothesis on M_1; value ℓ comes from forcing (set difference).

## Verification
- Exchange holds 134/134. The moved cell is the same cell in 240 steps, the φ-shifted cell in 208 steps, and some other cell in 0 steps. So the untranslated control fires.
- Chain invariant checked directly 582/582.
- SW: 987/987. The row variant gives 1426/1426 and is unused.
- MC: 1558/1558 (c2).
- Trustcheck is valid. A planted "G1 proved" is flagged on exactly MC and SW.

## Gaps
(MC) and (SW) are computed, not proved. The brief's premise for Task 1 ("local diamond, so finite case analysis") is **false**. A B/C interaction is a rerouting window of 2–7 hexagon steps, with 17 type-word classes on just 40 cells. Every B/C detour rejoins (779/779). No B/B detour rejoins (0/440), and that is the pre-registered control firing. A proof needs a uniform *crossing lemma*. Engine caveat: the intermediate composites lie outside Purbhoo's standard order. Both chain endpoints lie inside it.
