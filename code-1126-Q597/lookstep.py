import pickle, sys
from collections import Counter
def tris(S): return [t for t in S if t[0]=='tri']
def side(t):
    pts=list(t[1]); 
    if len(set(p[2:] for p in pts))==1: return 0, frozenset(p[:2] for p in pts), pts[0][2:]
    return 1, frozenset(p[2:] for p in pts), pts[0][:2]
def delta(s):
    """readout change of one step: {(side,tri): displacement}"""
    out={}
    for t in tris(s[1]):
        sd,k,w=side(t); out[(sd,k)]=[w,None]
    for t in tris(s[2]):
        sd,k,w=side(t); out.setdefault((sd,k),[None,None])[1]=w
    return {k:(tuple(b-a for a,b in zip(v[0],v[1])) if None not in v else v) for k,v in out.items()}
def fmt(d): return ' '.join('%d%s:%s'%(k[0],sorted(k[1])[0],v) for k,v in sorted(d.items(),key=str))
if __name__=='__main__':
    recs=pickle.load(open('/tmp/step0_5_%s.pkl'%sys.argv[1],'rb'))
    sys.path.insert(0,'.')
    for case,L,info in recs[:int(sys.argv[2])]:
        print('=== case',case,'L=',['x','y'][L])
        for r,(T,A,Ar,p,la,lb) in enumerate(info):
            print(' role',['x','y'][r],'T',T,'p',p,'mid',la,lb, 'L' if r==L else 'E')
            for nm,J in (('A ',A),('A\'',Ar)):
                for i in range(max(0,p-1),min(len(J),p+la+1)):
                    print('   ',nm,i,J[i][0],fmt(delta(J[i])))
