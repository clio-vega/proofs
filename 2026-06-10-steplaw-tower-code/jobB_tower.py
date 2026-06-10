"""
JOB B — map the Step-2 deflation / depth tower.

For each tie (|J*|>=2), the leading pi^mu coefficient cancels.  We chart the FULL
deflation: process terms by ascending val-level, record v_pi of the cumulative
partial sum at each level, and which indices activate at each level.  This is
carry-honest (we add exact Gaussian integers, never guess mod-2 digits).

  term_j = C(m,j) i^j (1+i)^j M_j,   v_pi(term_j) = val(j).
  level v_0=mu < v_1 < ...  ;  A_t = {j: val(j) <= v_t}  (cumulative active set).
  P_t = sum_{j in A_t} term_j  ;  trace = (v_t, newly-activated j, v_pi(P_t)).

  d = v_pi(G) - mu   (the survivor depth; +inf only for (2,2)).
  d1 = v_pi( sum_{j in J*} term_j ) - mu   (the within-J* / first-level depth).

Outputs results/step2-tower.csv and correlation verdicts.
Run:  python3 jobB_tower.py [MMAX]   (default 12)
"""
import sys, math, csv, os
from collections import defaultdict, Counter
from job1_tie_census import partitions, M_vector, v2, G_gaussian, vpi
from job2_mechanism import val_pieces, jstar_of, box_generators
from core_quotient import core_quotient

RESULTS = "/home/clio/projects/code/results"
os.makedirs(RESULTS, exist_ok=True)

def term_gaussian(m, j, Mj):
    """C(m,j) i^j (1+i)^j M_j as (re,im) exact integers."""
    coeff = math.comb(m, j) * Mj
    ij = [(1,0),(0,1),(-1,0),(0,-1)][j % 4]
    pr, pi = 1, 0
    for _ in range(j):
        pr, pi = pr - pi, pr + pi
    ar, ai = ij
    tr = ar*pr - ai*pi
    ti = ar*pi + ai*pr
    return coeff*tr, coeff*ti

def add(a, b):
    return (a[0]+b[0], a[1]+b[1])

def is_box(indices):
    """Is the index set an affine 2-adic box j0+{subset sums of distinct 2^a}?"""
    J = sorted(indices)
    if len(J) < 2:
        return (len(J)==1)  # singleton trivially a (0-dim) box
    if len(J) & (len(J)-1):  # size not a power of 2
        return False
    j0, S = box_generators(J)
    return S is not None

def deflation_trace(lam, m):
    Ms = M_vector(lam, m)
    # val per index (None if M_j=0)
    val = {}
    for j in range(m+1):
        if Ms[j] != 0:
            val[j] = j + 2*v2(math.comb(m,j)) + 2*v2(Ms[j])
    mu = min(val.values())
    J = sorted(j for j in val if val[j]==mu)
    # levels
    levels = sorted(set(val.values()))
    by_level = defaultdict(list)
    for j,v in val.items():
        by_level[v].append(j)
    # cumulative partial sums
    P = (0,0)
    trace = []   # (level, [new indices], vpi(P_cumulative))
    for v in levels:
        newj = sorted(by_level[v])
        for j in newj:
            P = add(P, term_gaussian(m, j, Ms[j]))
        trace.append((v, newj, vpi(*P)))
    # whole G
    reG, imG = G_gaussian(lam, m, Ms)
    vpG = vpi(reG, imG)
    d = None if vpG is None else vpG - mu
    # within-J* sum depth
    PJ = (0,0)
    for j in J:
        PJ = add(PJ, term_gaussian(m, j, Ms[j]))
    vpJ = vpi(*PJ)
    d1 = None if vpJ is None else vpJ - mu
    return dict(Ms=Ms, val=val, mu=mu, J=J, levels=levels, by_level=by_level,
                trace=trace, vpG=vpG, d=d, d1=d1, vpJ=vpJ)

def surviving_index_set(lam, m, info):
    """Indices j (val(j)<=vpG) whose own pi-adic digit at position vpG is 1.
    Carry-honest 'who is present at the surviving order'."""
    vpG = info['vpG']
    if vpG is None:
        return None
    Ms = info['Ms']; val = info['val']
    out = []
    for j,v in val.items():
        if v > vpG:
            continue
        re,im = term_gaussian(m, j, Ms[j])
        # pi-adic digit of term_j at position vpG = digit of unit at (vpG - v)
        # divide term by pi^v -> unit; then read bit at (vpG - v)
        for _ in range(v):
            # divide by (1+i): exact since v_pi(term)=v
            nre = (re+im)//2; nim=(im-re)//2; re,im=nre,nim
        t = vpG - v
        for _ in range(t):
            nre = (re+im)//2; nim=(im-re)//2; re,im=nre,nim
        if (re+im) % 2 == 1:
            out.append(j)
    return sorted(out)

def main(MMAX):
    rows=[]
    d_by_4core=defaultdict(set)
    d_by_4quot=defaultdict(set)
    d_by_size=defaultdict(Counter)
    survbox=0; survtot=0
    ddist=Counter()
    d1dist=Counter()
    within=0; needsnext=0
    for m in range(1,MMAX+1):
        for lam in partitions(2*m):
            info=deflation_trace(lam,m)
            if len(info['J'])<2:
                continue
            core,quot=core_quotient(lam,4)
            d=info['d']; d1=info['d1']
            surv=surviving_index_set(lam,m,info)
            ddist['inf' if d is None else d]+=1
            d1dist['inf' if d1 is None else d1]+=1
            if d is not None:
                d_by_4core[tuple(core)].add(d)
                d_by_4quot[tuple(tuple(q) for q in quot)].add(d)
                d_by_size[len(info['J'])][d]+=1
            # within-J* vs needs-next (does d==d1?)
            if d is not None and d1 is not None and d==d1:
                within+=1
            elif d is not None:
                needsnext+=1
            # surviving-order index set a box?
            if surv is not None:
                survtot+=1
                if is_box(surv):
                    survbox+=1
            rows.append(dict(
                m=m, lam='|'.join(map(str,lam)),
                Jstar='|'.join(map(str,info['J'])), Jsize=len(info['J']),
                mu=info['mu'],
                vpi=('inf' if info['vpG'] is None else info['vpG']),
                d=('inf' if d is None else d),
                d1=('inf' if d1 is None else d1),
                surv_set=('-' if surv is None else '|'.join(map(str,surv))),
                surv_size=(0 if surv is None else len(surv)),
                surv_is_box=('-' if surv is None else int(is_box(surv))),
                core4='|'.join(map(str,core)),
                quot4=';'.join('|'.join(map(str,q)) if q else '' for q in quot),
                parity=info['J'][0]%2,
                v2f=v2(info['Ms'][0]),
                nlevels=len(info['levels']),
            ))
    with open(f"{RESULTS}/step2-tower.csv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote results/step2-tower.csv ({len(rows)} ties, m<={MMAX})")
    print(f"\nd = v_pi(G)-mu distribution: {dict(sorted(ddist.items(),key=lambda kv:str(kv[0])))}")
    print(f"d1 = within-J* depth distribution: {dict(sorted(d1dist.items(),key=lambda kv:str(kv[0])))}")
    print(f"\nsurvivor WITHIN J* (d==d1): {within}")
    print(f"survivor NEEDS next-val levels (d>d1): {needsnext}")
    print(f"surviving-order index set is an affine 2-adic box: {survbox}/{survtot}")

    # correlation: is d a function of 4-core alone? 4-quotient alone?
    sc=sum(1 for v in d_by_4core.values() if len(v)==1)
    sq=sum(1 for v in d_by_4quot.values() if len(v)==1)
    print(f"\nd determined by 4-core alone:     {sc}/{len(d_by_4core)} single-valued")
    print(f"d determined by 4-quotient alone: {sq}/{len(d_by_4quot)} single-valued")
    print(f"\nd by |J*|:")
    for sz in sorted(d_by_size):
        print(f"  |J*|={sz}: {dict(sorted(d_by_size[sz].items()))}")

if __name__=="__main__":
    MMAX=int(sys.argv[1]) if len(sys.argv)>1 else 12
    main(MMAX)
