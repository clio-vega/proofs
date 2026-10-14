"""Task 2 / Q457: degrees of c_{(1^n)} at k>=2, far beyond engine M1's reach."""
import sys, time
sys.path.insert(0,'/home/clio/projects/proofs/code-1014-q458')
from censusk import c_1n

def run(k, ell, nmax):
    M = k+ell+1; N = 2*M
    print(f"### k={k} ell={ell}  M={M} N={N}  (h=2k={2*k}, w=2ell+2={2*ell+2})", flush=True)
    for n in range(2, nmax+1, 2):
        t0=time.time()
        d = c_1n(n,k,ell)
        if not d:
            print(f"  n={n:3d}: c=0  ({time.time()-t0:.1f}s)", flush=True); continue
        deg=max(d); lo=min(d)
        # k=1 prediction, applied verbatim to k>=2 (expected to FAIL -> negative control)
        ys=[y for y in range(-(n+4),n+5) if -n<=N*y<=n+2]
        pred=max(M*y*y-y for y in ys)
        print(f"  n={n:3d}: deg={deg:4d} lowdeg={lo:3d} nterms={len(d)} "
              f"k1law={pred:4d} {'MATCH' if deg==pred else 'differ'}  n^2/(4M)={n*n/(4*M):.1f}"
              f"  ({time.time()-t0:.1f}s)", flush=True)

if __name__=='__main__':
    for k,ell,nmax in [(2,0,28),(2,1,28),(3,0,24)]:
        run(k,ell,nmax)
