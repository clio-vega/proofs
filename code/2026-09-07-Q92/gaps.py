import sys, itertools, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B, engineC as C
t=sp.Symbol('t')

print("=== GAP 1: uniform witness  lambda=empty, b=-1, c=-f  for all e,f>=2 ===")
print("   predicted two-bead matrix element:  t^(f-1) (t - t^-1) = t^f - t^(f-2)")
ok=bad=0
for e in range(2,8):
    for f in range(2,8):
        if e==f: continue
        M=B.maya(())
        b,c=-1,-f
        Mp = M-{b,c}|{b+e,c+f}
        mu=B.unmaya(Mp)
        pred = t**f - t**(f-2)
        act  = A.commutator((),e,f).get(mu,sp.Integer(0))
        # confirm it really is a two-bead target and sites are distinct
        distinct = len({b,c,b+e,c+f})==4
        twobead  = len(M^Mp)==4
        good = distinct and twobead and sp.simplify(act-pred)==0
        ok+=good; bad+= (not good)
        if not good: print("  FAIL",e,f,mu,"pred",pred,"act",act,distinct,twobead)
print(f"   {ok} confirmed, {bad} failed   (2<=e,f<=7, e!=f)")

print("\n=== GAP 2: factorization of the 3-bead sector of [R_e,[R_f,R_g]] ===")
print("   predicted:  t^W (t^k32 - t^-k32)(t^(k21+k31) - t^-(k21+k31))  summed over assignments")
def kap(bi,ei,bj,ej):
    return (1 if bj<bi+ei<bj+ej else 0)-(1 if bj<bi<bj+ej else 0)
def triple_comm(lam,e,f,g):
    v={lam:sp.Integer(1)}
    def op(vec,n): return A.apply_op(vec,n)
    def sub(a,b):
        o={}
        for k in set(a)|set(b):
            x=sp.expand(a.get(k,0)-b.get(k,0))
            if x!=0: o[k]=x
        return o
    inner_on=lambda w: sub(op(op(w,4 if False else g),f), op(op(w,f),g))  # [R_f,R_g]w
    return sub(op(inner_on(v),e), inner_on(op(v,e)))
def predict3(lam,e,f,g,mu):
    M=B.maya(lam); Mp=B.maya(mu)
    R=sorted(M-Mp); Ad=sorted(Mp-M)
    if len(R)!=3: return None
    tot=sp.Integer(0)
    for perm in itertools.permutations(Ad):
        # bead R[i] moves by disp[i]; need multiset of disps == {e,f,g} with the
        # assignment b1->e, b2->f, b3->g
        for assign in itertools.permutations([(e,0),(f,1),(g,2)]):
            pass
        break
    # enumerate assignments directly: choose which removed bead takes e, f, g
    for pr in itertools.permutations(R):
        b1,b2,b3=pr
        if sorted([b1+e,b2+f,b3+g])!=Ad: continue
        # legality of every intermediate is implied: all six sites distinct
        if len({b1,b2,b3,b1+e,b2+f,b3+g})!=6: continue
        P=(C.cnt(M,b1,b1+e)+C.cnt(M,b2,b2+f)+C.cnt(M,b3,b3+g))
        k32=kap(b3,g,b2,f); k21=kap(b2,f,b1,e); k31=kap(b3,g,b1,e)
        K=k21+k31
        tot+=sp.expand(t**P*(t**k32-t**(-k32))*(t**K-t**(-K)))
    return sp.expand(sp.cancel(tot))
def nm(lam,mu): return len(B.maya(lam)^B.maya(mu))//2
ok=bad=0; shown=0
for (e,f,g) in [(2,3,4),(2,4,3),(3,2,4),(2,3,5),(3,4,5),(2,5,3)]:
    for lam in [(),(1,),(2,),(1,1)]:
        out=triple_comm(lam,e,f,g)
        for mu,co in out.items():
            if nm(lam,mu)!=3: continue
            pr=predict3(lam,e,f,g,mu)
            if pr is None: continue
            if sp.simplify(co-pr)==0: ok+=1
            else:
                bad+=1
                if shown<5: print("  FAIL",(e,f,g),lam,mu,"actual",co,"pred",pr); shown+=1
print(f"   {ok} agree, {bad} disagree   (3-bead matrix elements)")
