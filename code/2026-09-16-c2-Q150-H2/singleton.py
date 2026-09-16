import sys
from twobead import words, relations_2Bstar, relations_2Bstarstar, support_ok

def borders(w):
    n=len(w); return [p for p in range(n) if w[:p]==w[n-p:]]
def crit_star(w):
    n=len(w); return all(w[p] != w[n-1-p] for p in borders(w))
def crit_starstar(w):
    n=len(w); d=n+1
    return all(w[dd-1] != w[d-dd-1] for dd in range(1,d))

for d in range(2, 10):
    ws = words(d-1)
    rs  = list(relations_2Bstar(d))
    rss = list(relations_2Bstarstar(d))
    ok_s  = set(w for w in ws if support_ok(frozenset([w]), rs))
    ok_ss = set(w for w in ws if support_ok(frozenset([w]), rss))
    p_s   = set(w for w in ws if crit_star(w))
    p_ss  = set(w for w in ws if crit_starstar(w))
    print('d=%2d  (2B*): actual %5d pred %5d match=%-5s |  (2B**): actual %5d pred %5d match=%s'
          % (d, len(ok_s), len(p_s), ok_s==p_s, len(ok_ss), len(p_ss), ok_ss==p_ss))
    sys.stdout.flush()
