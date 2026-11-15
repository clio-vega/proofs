"""Exchange test (Q585 G1).  M = census mosaic; (ii) = flock sgn=+1 on A, migrate N (to B); (iii) = flock sgn=-1, migrate W (to C).
M_1 = M after (ii) has sent ONLY its first block = the top row of alpha (max flock value, processed first).
(Exchange) landings of (iii) on M_1 (as raw C-nest positions) == first m - alpha_top landings of (iii) on M.
Also records the SOURCE cells (raw A-nest positions) of each (iii) step on M and on M_1, to test the translation conjecture."""
import sys, copy, pickle
sys.path.insert(0,'/home/clio/projects/proofs/code-1114c2-Q590')
from common import *
from readouts import unique_flock
from collections import Counter
def run_iii(b):
    """(iii) on board b (modified in place): returns list of (source A-pos, landing C-pos)."""
    u, w = b.nest[0]['u'], b.nest[0]['w']; E_ = tuple(-a for a in u); N_ = tuple(-a for a in w)
    fl = unique_flock(b, 0, E_, N_)
    dvec = tuple(-a for a in E_)          # 'W' relative to (-E,-N) = +u = E_A
    key = {t: (v, box_of(b, 0, b.tiles[t][1], E_, N_)[0]) for t, v in fl.items()}
    order = sorted(fl, key=lambda t: key[t])
    inv = {t: ij for ij, t in b.nest[0]['pos'].items()}; out = []
    for t in order:
        del b.nest[0]['pos'][inv[t]]; R, X, ij, s = b.journey(t, 0, unit_dir(dvec))
        assert X == 2; out.append((inv[t], ij, fl[t]))
    return out
def run_ii_prefix(b, top_only=True):
    u, w = b.nest[0]['u'], b.nest[0]['w']; E_, N_ = u, w
    fl = unique_flock(b, 0, E_, N_)
    key = {t: (v, box_of(b, 0, b.tiles[t][1], E_, N_)[0]) for t, v in fl.items()}
    order = sorted(fl, key=lambda t: key[t], reverse=True)
    vmax = max(fl.values()); inv = {t: ij for ij, t in b.nest[0]['pos'].items()}; out = []
    for t in order:
        if top_only and fl[t] != vmax: break
        del b.nest[0]['pos'][inv[t]]; R, X, ij, s = b.journey(t, 0, unit_dir(N_))
        assert X == 1; out.append((inv[t], ij, fl[t]))
    return out
if __name__ == '__main__':
    st = Counter(); recs = []; first_bad = None
    for key, b0 in census():
        A0 = sorted(b0.nest[0]['pos'])
        if len(set(j for i, j in A0)) < 2 and len(set(i for i, j in A0)) < 2: pass
        M_iii = run_iii(copy.deepcopy(b0))
        b1 = copy.deepcopy(b0); top = run_ii_prefix(b1)
        if not b1.nest[0]['pos']: st['alpha is one row (vacuous)'] += 1; continue
        M1_iii = run_iii(b1)
        k = len(M1_iii); st['n'] += 1
        landM = [l for s, l, v in M_iii][:k]; landM1 = [l for s, l, v in M1_iii]
        ok = landM == landM1; st['EXCHANGE holds'] += ok
        srcM = [s for s, l, v in M_iii][:k]; srcM1 = [s for s, l, v in M1_iii]
        st['sources identical (untranslated)'] += srcM == srcM1
        recs.append(dict(key=key, A=A0, top=[s for s, l, v in top], M_iii=M_iii, M1_iii=M1_iii, ok=ok))
        if not ok and first_bad is None:
            first_bad = (key, A0, [s for s, l, v in top], M_iii, M1_iii)
    print(dict(st)); print('first mismatch:', first_bad)
    pickle.dump(recs, open('/tmp/c2-exchange.pkl', 'wb'))
