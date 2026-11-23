"""Lesion of Lemma U inside the local classification: e(s) on the survivors of ccloc with LESION noV (drop the re-run half
of Theorem fc, i.e. uniqueness on the re-run's board) and noS2 (drop the choice rule J2 and UQ, keep V)."""
import pickle
exec(open('/home/clio/projects/proofs/code-1123-Q597/commonpot.py').read().split("st = Counter()")[0])
for tag in ('-noV', '-noS2'):
    st = Counter()
    for s in pickle.load(open('/home/clio/projects/proofs/code-1123-Q597/ccsurv%s.pkl' % tag, 'rb')):
        T = s['T']; To = OTHER[T]; Rc = s['seg'][s['c0']][0]; R1 = step(T, s['RA'], s['t'])[1]
        e = sub(Pi(Rc, To), Pi(R1, T)); rd = sorted(u[0] for u in s['seg'][s['c0']][2] if u in set(s['S']))
        st[(tag, T, 'e<=0', all(v <= 0 for v in e), 'e!=0', any(e), 'c0 reads', tuple(rd), 'J c0', (cls(To, Rc), s['seg'][s['c0']][1]))] += 1
    for k in sorted(st, key=str): print(st[k], k)
