# Q590 — Order Reversal vs switching growth diagram (PROVE 2026-11-14 c2)
**Statement tested.** Is Order Reversal (c_k = ι(b_{m+1−k}), equivalent to Purbhoo Prop 5.1 (ii)=(iii)) one switching growth diagram, as in Lemma D?
**Result.**
- The literal form is ill-posed: the two named chains are equal.
- The well-posed form, single-migration growth over J(α), is **refuted**: 233/1100 ideals for B, 129/1100 for C, with an explicit witness.
- Step 0: the control is the target. The identity is not bookkeeping: step strings match on 86/606 journeys, and move counts differ, 4262 vs 3612.
- Computed: **mixed commutation**. A→B and A→C single migrations commute 1558/1558, including 262 interacting cases. B/B pairs commute in 0 of their 440 interacting cases.
- Commutation does not imply Order Reversal: hypothesis K fails 342/740.
- The gap is now located at "Exchange", in the induction on ℓ(α).
- The coplactic lemma locator is van Leeuwen math/9908099 §3.3, quoted verbatim.
**Gaps.** G1 is not proved. Mixed commutation is not proved. Exchange is open.
Full write-up: `2026-11-14c2-Q590-order-reversal-switching.tex` (commit 3fb80b68369e; code @ 64131e56a304).
