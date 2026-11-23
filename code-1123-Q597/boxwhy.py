"""Why the failure box is empty: per survivor pair, the lambda and mu slacks.  Also the non-strict box (R^y_{p_y} in place
of R^y_{p_y+1}, the Share case) as a control: it must contain the D of every Share-consistent pair."""
import os
exec(open('/home/clio/projects/proofs/code-1123-Q597/pairloc.py').read().split("\nZ = (0, 0, 0, 0)\n")[0])
st = Counter()
Z = (0, 0, 0, 0)
def sub(a, b): return tuple(p - q for p, q in zip(a, b))
for sx in SB:
    x = role(sx, 'B', Z)
    for sy in SC:
        y0 = role(sy, 'C', Z)
        Ryc, Rxc = x['seg'][0][0], y0['seg'][0][0]
        # lambda: a_x = lam(Ryc)-lam(RAx), a_y = lam(Rxc)-lam(RAy0) (frame-free); delta = lam(R1)-lam(RA)
        ax, ay = sub(lam(Ryc), lam(x['RA'])), sub(lam(Rxc), lam(y0['RA']))
        dx, dy = sub(lam(x['R1']), lam(x['RA'])), sub(lam(y0['R1']), lam(y0['RA']))
        slack = tuple(a + b - c - d for a, b, c, d in zip(ax, ay, dx, dy))     # failure needs slack >= 0 componentwise
        nonstrict = tuple(a + b for a, b in zip(ax, ay))                           # Share-type: needs >= 0
        st[('lam: a_x', ax, 'a_y', ay, 'd_x', dx, 'd_y', dy, 'fail-slack', slack, 'share-slack', nonstrict)] += 1
        st[('lam fail-slack min', min(slack))] += 1
        # share control: D forced, then lambda equalities must hold if both segment starts are the contact steps
        D = sub(base(Ryc), base(y0['RA']))
        if tr(y0['RA'], D) == Ryc and Ryc and x['seg'][0][1] == sy['t']:
            y = role(sy, 'C', D)
            consistent = y['seg'][0][0] == x['RA'] and y['seg'][0][1] == sx['t']
            st[('SHARE-consistent pair', consistent, 'lam share-slack', nonstrict)] += 1
for k in sorted(st, key=str): print(st[k], k)
