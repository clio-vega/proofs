"""Verify STEP 1 of the proof: on Gamma, the ONLY q that can witness a symmetric
exchange for p is q = n+p (p in the x-block) resp. q = p-n (p in the y-block).
This tests the four-case analysis directly, not just its conclusion."""
import itertools, random
random.seed(7)
checked = succ = 0
bad_case = 0
for trial in range(20000):
    n = random.choice([1,2,3]); m = random.choice([1,2,3])
    pts = list(itertools.product(range(m+1), repeat=n))
    A = set(random.sample(pts, random.randint(1, min(len(pts), 5))))
    S = {tuple(a) + tuple(m-ai for ai in a) for a in A}
    N = 2*n
    for u in S:
        for v in S:
            if u == v: continue
            for p in range(N):
                if u[p] <= v[p]: continue
                checked += 1
                for q in range(N):
                    if u[q] >= v[q]: continue
                    a2 = list(u); a2[p]-=1; a2[q]+=1
                    b2 = list(v); b2[q]-=1; b2[p]+=1
                    if tuple(a2) in S and tuple(b2) in S:
                        succ += 1
                        partner = p+n if p < n else p-n
                        if q != partner:
                            bad_case += 1
                            if bad_case <= 3:
                                print("  UNEXPECTED WITNESS n=%d m=%d p=%d q=%d u=%s v=%s"%(n,m,p,q,u,v))
print("(p,u,v) triples with u_p>v_p examined : %d" % checked)
print("  of these, SOME q witnesses exchange : %d" % succ)
print("  witness q != the forced partner     : %d   <-- must be 0" % bad_case)

# and: verify that on Gamma, u-e_p+e_q lands back on Gamma ONLY for q=partner
off = 0; on = 0
for trial in range(3000):
    n = random.choice([1,2,3]); m = random.choice([1,2,3])
    a = tuple(random.randint(0,m) for _ in range(n))
    u = a + tuple(m-ai for ai in a); N = 2*n
    for p in range(N):
        for q in range(N):
            if p == q: continue
            w = list(u); w[p]-=1; w[q]+=1
            onG = all(w[i]+w[n+i] == m for i in range(n))
            partner = p+n if p < n else p-n
            if onG: on += 1
            if onG != (q == partner): off += 1
print("\nlattice check: u-e_p+e_q lands on Gamma  <=>  q = partner(p)")
print("  landings on Gamma : %d   violations of the equivalence : %d  <-- must be 0" % (on, off))
