"""KEY LEMMA: if beta = c_+ for c concave on Z with bounded positive support,
then gamma(y) = sum_{r>=y} beta(r) is log-concave.

Instrument discipline: this script must (i) CONFIRM the lemma on truncated-concave
beta, and (ii) REFUSE it -- i.e. find failures -- when beta is merely PF2/log-concave
but not truncated-concave.  An engine that reports 0 failures on both is blind."""
import random, itertools

def gamma_of(beta, lo):
    """beta a list starting at index lo. Return dict y->tail sum."""
    g={}; acc=0
    for k in range(len(beta)-1,-1,-1):
        acc+=beta[k]; g[lo+k]=acc
    return g

def tail_logconcave(beta, lo):
    g=gamma_of(beta,lo)
    get=lambda y: g.get(y, 0 if y>lo+len(beta)-1 else sum(beta))
    bad=[]
    for y in range(lo-1, lo+len(beta)+2):
        G=get(y); Gm=get(y-1); Gp=get(y+1)
        if G*G < Gm*Gp: bad.append((y,G,Gm,Gp))
    return bad

def is_concave_seq(c):
    return all(2*c[i]>=c[i-1]+c[i+1] for i in range(1,len(c)-1))
def is_lc_pos(b):
    # log-concave with interval support
    nz=[i for i,v in enumerate(b) if v>0]
    if nz and nz[-1]-nz[0]+1!=len(nz): return False
    return all(b[i]**2>=b[i-1]*b[i+1] for i in range(1,len(b)-1))

random.seed(20261001)
print("=== (i) CONFIRM: beta = positive part of a concave integer sequence ===")
fails=0; trials=0
for _ in range(300000):
    Lw=random.randint(2,9)
    # build a concave integer sequence by choosing non-increasing increments
    start=random.randint(-6,8)
    incs=sorted([random.randint(-5,5) for _ in range(Lw-1)], reverse=True)
    c=[start]
    for d in incs: c.append(c[-1]+d)
    beta=[max(v,0) for v in c]
    if sum(beta)==0: continue
    trials+=1
    bad=tail_logconcave(beta,0)
    if bad:
        fails+=1
        if fails<=3: print("   FAIL",c,beta,bad)
print(f"   trials={trials} failures={fails}")

print()
print("=== (ii) REFUSE: beta merely log-concave (PF2) but NOT truncated-concave ===")
fails=0; trials=0; shown=0
for _ in range(400000):
    Lw=random.randint(3,7)
    beta=[random.randint(0,9) for _ in range(Lw)]
    if sum(beta)==0: continue
    if not is_lc_pos(beta): continue
    nz=[i for i,v in enumerate(beta) if v>0]
    sup=beta[nz[0]:nz[-1]+1]
    if is_concave_seq(sup): continue      # exclude the truncated-concave ones
    trials+=1
    bad=tail_logconcave(beta,0)
    if bad:
        fails+=1
        if shown<4: print("   refused (as it must):",beta,bad[:2]); shown+=1
print(f"   PF2-but-not-concave trials={trials} failures={fails}")
if trials and fails==0: print("   *** WARNING: instrument may be blind -- no refusal found ***")

print()
print("=== (iii) exhaustive small check: all concave c with values in [-3,6], length<=6 ===")
fails=0; trials=0
for Lw in range(1,7):
    for c in itertools.product(range(-3,7), repeat=Lw):
        if not is_concave_seq(list(c)): continue
        beta=[max(v,0) for v in c]
        if sum(beta)==0: continue
        trials+=1
        if tail_logconcave(beta,0):
            fails+=1
            if fails<=3: print("   FAIL",c)
print(f"   exhaustive trials={trials} failures={fails}")
