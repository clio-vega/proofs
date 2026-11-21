"""PROVE 11-21 c2: the lesions of Step 0 recomputed at J_o-STEP level from the typed words in sheetvac5.out
(the first in-run lesion used the raw step type, which is 'a'/'b', so '#a1' was identically 0: instrument error).
Statistic per starred (read) J_o step: sheet = #b before mod 2 / excursion form (not a class-1 a-step);
lesion L1 = #a1 before mod 2; lesion L2 = j mod 2; lesion L3 = #a0 before mod 2.  usage: python3 lesion_words.py"""
import ast
from collections import Counter
st = Counter()
for line in open('sheetvac5.out'):
    if 'WORD Jo' not in line: continue
    n, rest = line.split(' ', 1); n = int(n); t = ast.literal_eval(rest); tag, w = t[0], t[2].split()
    typ = [x.rstrip('*') for x in w]
    for j, x in enumerate(w):
        if not x.endswith('*'): continue
        c, ty = int(typ[j][0]), typ[j][1:]
        sheet = sum(1 for z in typ[:j] if z.endswith('b')) % 2
        assert sheet == c, (w, j)           # precondition: class field == #b parity
        L1 = sum(1 for z in typ[:j] if z.endswith('a1')) % 2
        L2 = j % 2
        L3 = sum(1 for z in typ[:j] if z.endswith('a0')) % 2
        st[(tag, 'sheet law (not class-1 a-step)', not (c == 1 and ty != 'b'))] += n
        st[(tag, 'L1 #a1 parity', L1, 'differs from sheet_pre', L1 != sheet)] += n
        st[(tag, 'L2 j parity', L2, 'differs from sheet_pre', L2 != sheet)] += n
        st[(tag, 'L3 #a0 parity', L3, 'differs from sheet_pre', L3 != sheet)] += n
for k in sorted(st, key=str): print(st[k], k)
