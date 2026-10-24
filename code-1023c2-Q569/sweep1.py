from puzzles import *
from kalg import *
import itertools, sys

KTW   = ['zero_A','zero_B','one_A','one_B','rh_014','rh_230','rh_125']
SET_K   = KTW + ['K']
SET_TRI = KTW + ['tri']
SET_HEX = KTW + ['hex']
SET_HEXR= KTW + ['hexR']

def parts_in(rows, cols):
    out=[]
    def rec(pref, prev, r):
        if r==0: out.append(tuple(pref)); return
        for v in range(min(prev,cols), -1, -1):
            rec(pref+[v], v, r-1)
    rec([], cols, rows)
    return [tuple(x for x in p if x>0) for p in out]

def pad(p,k): return list(p)+[0]*(k-len(p))

n,k = int(sys.argv[1]), int(sys.argv[2])
which = sys.argv[3]
NK = n-k
print("### n=%d k=%d  ambient = %d rows x %d cols" % (n,k,k,NK))

if which=='c':      # Vakil \K  -> c   (ambient k x (n-k)), |nu|>=|la|+|mu|
    P = parts_in(k, NK); tiles=SET_K; mx=k*NK
    rows=0; bad=0
    for la,mu in itertools.product(P,P):
        cc = c_coeffs(la,mu,max(4,k+len(mu)+1),maxdeg=mx)
        for nu in P:
            if sum(nu) < sum(la)+sum(mu): continue
            exp = cc.get(tuple(x for x in nu if x>0),0)
            got = count_puzzles(n,pad(la,k),pad(mu,k),pad(nu,k),k,tiles)
            rows+=1
            if exp!=got:
                bad+=1; print("  MISMATCH la=%s mu=%s nu=%s  alg=%d puzzle=%d"%(la,mu,nu,exp,got))
    print("  \\K vs c      : %d/%d agree"%(rows-bad,rows))

if which=='ct':     # WZ \tri -> c~
    P = parts_in(k, NK); tiles=SET_TRI; mx=k*NK
    rows=0; bad=0
    for la,mu in itertools.product(P,P):
        cc = ct_coeffs(la,mu,max(4,k+len(mu)+1),maxdeg=mx)
        for nu in P:
            if sum(nu) < sum(la)+sum(mu): continue
            exp = cc.get(tuple(x for x in nu if x>0),0)
            got = count_puzzles(n,pad(la,k),pad(mu,k),pad(nu,k),k,tiles)
            rows+=1
            if exp!=got:
                bad+=1; print("  MISMATCH la=%s mu=%s nu=%s  alg=%d puzzle=%d"%(la,mu,nu,exp,got))
    print("  \\tri vs c~   : %d/%d agree"%(rows-bad,rows))

if which=='d':      # \hex -> d   (ambient k x (n-k-1)), |nu|<=|la|+|mu|
    P = parts_in(k, NK-1); tiles=SET_HEX
    rows=0; bad=0; dc={}
    for nu in P:
        dc[nu] = d_coeffs(nu, max(2,k), max(2,k), maxdeg=2*max(1,sum(nu)))
    for la,mu in itertools.product(P,P):
        for nu in P:
            if sum(nu) > sum(la)+sum(mu): continue
            exp = dc[nu].get((tuple(x for x in la if x>0),tuple(x for x in mu if x>0)),0)
            got = count_puzzles(n,pad(la,k),pad(mu,k),pad(nu,k),k,tiles)
            rows+=1
            if exp!=got:
                bad+=1; print("  MISMATCH la=%s mu=%s nu=%s  alg=%d puzzle=%d"%(la,mu,nu,exp,got))
    print("  \\hex vs d    : %d/%d agree   (domain: la,mu,nu in %d rows x %d cols)"%(rows-bad,rows,k,NK-1))

if which=='dt':     # \hexR -> d~ = d^{nu'}_{la' mu'}  (ambient (k-1) x (n-k))
    P = parts_in(k-1, NK); tiles=SET_HEXR
    rows=0; bad=0; dc={}
    for nu in P:
        nt=tuple(transpose(nu))
        dc[nu] = d_coeffs(nt, max(2,NK), max(2,NK), maxdeg=2*max(1,sum(nu)))
    for la,mu in itertools.product(P,P):
        for nu in P:
            if sum(nu) > sum(la)+sum(mu): continue
            key=(tuple(transpose(la)),tuple(transpose(mu)))
            exp = dc[nu].get(key,0)
            got = count_puzzles(n,pad(la,k),pad(mu,k),pad(nu,k),k,tiles)
            rows+=1
            if exp!=got:
                bad+=1; print("  MISMATCH la=%s mu=%s nu=%s  alg=%d puzzle=%d"%(la,mu,nu,exp,got))
    print("  \\hexR vs d~  : %d/%d agree   (domain: la,mu,nu in %d rows x %d cols)"%(rows-bad,rows,k-1,NK))
