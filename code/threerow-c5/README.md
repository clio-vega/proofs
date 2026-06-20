# three-row c=5 interior even-|J*| (PROVE 2026-06-19)

Closed form, structural decomposition, and the interior even-|J*| theorem for (a,b,5).
Engine reuses ../threerow-c3/mn.py (Murnaghan-Nakayama).

## Result
For box-interior (a,b,5) the interior minimiser set J* (only generator 2 fires, |J*|<=2):
- a even (b odd): j0=0, J*={0,2} iff (a,b)=(0,1) mod 4, else {0}.
- a odd  (b even): j0=3, J*={3,5} iff (a,b)=(3,0) mod 4, else {3}.

## Engine
Q_5(j) = (a-3)(b-4) H_5(j) - 10! C(j,10),  H_5 = sum_{k=0}^8 h_k C(j,k) (h_k saved in H5sym.pkl)
H_5(0) = (a+3)(a+4)(a+5)(a+6)(b+2)(b+3)(b+4)(b+5)
c=5 Number Lemma: 8 | H_5 (v2 H_5 >= 3), constant floor, sharp at (a,b,j)=(3,0,2).

## scripts
census.py census2.py          empirical J* distribution
closedform.py H5.py prop2.py  closed form Q_5, H_5 reconstruction, Prop2 Kummer (0 mism)
peel_deep.py                  peeling identity (0 mism) + g(j) exceptions {10,11,16,17,18,19}
exceptional_rig.py            deep exceptional indices: 2^{w_j}|C C H5 (0 fail K=9,11)
lowfinal2.py                  j=4..9: v2(Psi_j)>=kappa_j (0 fail K=8,10)
lowj123.py                    j=1,2,3 reductions + R-valuations + tie classes
endtoend.py                   full theorem vs MN: 561 shapes, 0 mism (m<=45)
