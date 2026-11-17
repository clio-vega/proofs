import pickle
from collections import Counter
S = pickle.load(open('rolesigs.pkl', 'rb'))
NAMES = ['U0=up(u)', 'U1=up(u-e1)', 'U2=up(u-e2)', 'D1=dn(u-e1)', 'D2=dn(u-e2)', 'D3=dn(u-e1-e2)']
sym = {0: {(0, 0): '0', (0, -1): 's', (1, -1): 'b', (-1, 0): 'c'}, 1: {(0, 0): '0', (0, 1): 's', (-1, 1): 'b', (1, 0): 'c'}}
for (side, role), ss in sorted(S.items()):
    if not ss: continue
    (sg,) = ss
    print('side', 'q' if side == 0 else 'p', role.ljust(6), '  '.join(f'{n}:({sym[side][x]},{sym[side][y]})' for n, (x, y) in zip(NAMES, sg)))
