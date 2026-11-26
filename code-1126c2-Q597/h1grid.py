"""(H1) in all cases, and non-stop of E^L after E^L_0, by the stop grid (casea tex Lemma grid, Lemma gridmono table).
Inputs per class: E_0 (Lemma loc2), E^L_0 and E^L_1 (Prop rerun: presence + Lemma U; E^L_1 needs E^L not to stop at R'),
J_E's known steps.  Order of the argument (no circularity):
  T_E = C:  [E^L non-stop at R'] from E_0 alone  ->  E^L_1 forced  ->  [J_E non-stop at R_E^1] from E^L_0 + E^L_1.
  T_E = B:  [J_E non-stop at R_E^1] from E^L_0 alone  ->  (H1)  ->  [E^L non-stop at R'] from E_0 + E_1 (+ E_m=(0,b) in case b, by F3).
Monotone: B non-decreasing, C non-increasing.  'X stops at R' => F(R) = c (nest T_E grid; same nest on M0 and M0.J_L)."""
dF = {'B': {(0,'a0'):(0,0),(0,'a1'):(0,1),(0,'b'):(1,1),(1,'a0'):(1,0),(1,'a1'):(0,1),(1,'b'):(0,0)},
      'C': {(0,'a0'):(-1,0),(0,'a1'):(0,0),(0,'b'):(0,0),(1,'a0'):(-1,0),(1,'a1'):(0,-1),(1,'b'):(-1,-1)}}
aT = {'B':'a1','C':'a0'}
def tot(T, steps): return tuple(sum(dF[T][s][i] for s in steps) for i in range(2))
def le(a, b): return all(x <= y for x, y in zip(a, b))
def excluded(T, stop_steps, other_steps):
    """journey X stops after stop_steps (so c - F0 = tot(stop)); journey Y from the same rhombus R_E took other_steps
    (a prefix of its steps from R_E) and ends at c: B needs tot(other) <= c - F0, C needs tot(other) >= c - F0.
    Returns True iff that is impossible."""
    s, o = tot(T, stop_steps), tot(T, other_steps)
    return not le(o, s) if T == 'B' else not le(s, o)
ok = True
for TL, E0, EL0 in (('B',(1,'b'),(1,'a1')),('B',(0,'a0'),(0,'b')),('C',(0,'a1'),(0,'b')),('C',(1,'b'),(1,'a0'))):
    TE = {'B':'C','C':'B'}[TL]; EL1 = (1, 'b')
    for case in 'ab':
        E1 = (0, aT[TL]) if case == 'b' else (0, 'b')
        if TE == 'C':
            a = excluded(TE, [EL0], [E0])            # E^L stops at R'?  J_E from R_E took E_0
            b = excluded(TE, [E0], [EL0, EL1])       # J_E stops at R_E^1?  E^L took EL0, EL1
        else:
            b = excluded(TE, [E0], [EL0])            # J_E stops at R_E^1?
            a = excluded(TE, [EL0], [E0, E1] + ([(0, 'b')] if case == 'b' else []))   # E^L stops at R'? uses E_1 (after H1) and, in case b, E_m=(0,b) by (F3)
        ok &= a and b
        print('T_E', TE, 'E0', E0, 'EL0', EL0, 'EL1', EL1, 'case', case, '| E^L cannot stop at R\':', a, '| J_E cannot stop at R_E^1 (H1):', b)
print('ALL EXCLUDED:', ok)
