"""Independent KT puzzle enumerator (labels 0,1,2:=10). Row r from bottom; up(r,j) j<N-r: bottom h[r][j], left a[r][j] (/), right b[r][j] (\\).
down(r,j) j<N-r-1: top h[r+1][j], left(\\) b[r][j], right(/) a[r][j+1].  NW boundary a[r][0] (r bottom->top), NE boundary b[r][N-r-1]."""
import sys; sys.path.insert(0,'/home/clio/git/puzzles')
DELTA={(0,0,0),(1,0,2),(2,1,0),(1,1,1),(0,2,1)}
NAB={(0,0,0),(0,2,1),(2,1,0),(1,0,2),(1,1,1)}
def puzzles(S,NW,NE,nab):
    N=len(S); out=[]
    def row(r,h,acc):
        w=N-r
        if w==0: out.append(acc); return
        a0=NW[r]
        def go(j,a,b,hn):
            # up(r,j): bottom h[j], left a[j]
            outs=[x for x in (0,1,2) if (h[j],a[j],x) in DELTA]
            for bj in outs:
                if j==w-1:
                    if bj!=NE[r]: continue
                    row(r+1,hn,acc+[(tuple(a),tuple(b+[bj]),tuple(hn))]); continue
                for top in (0,1,2):
                    for ar in (0,1,2):
                        if nab(top,bj,ar):
                            go(j+1,a+[ar],b+[bj],hn+[top])
        go(0,[a0],[],[])
    row(0,list(S),[])
    return out
readings={'R1 (top,left\\,right/)':lambda t,l,r:(t,l,r) in NAB,'R2 (top,right/,left\\)':lambda t,l,r:(t,r,l) in NAB}
