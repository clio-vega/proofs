"""strict mode: (i) positive control: A half crossing exactly ONE B half once (and the other B half
untouched) -- the one-crossing pictures of Rem 10.3 must exist; (ii) pigeonhole-reach statistic:
A half crosses one B half once and its rhombus chain meets the OTHER B half in a shared rhombus where
both strands go straight through the SAME edge class (a same-class collision = the pigeonhole step)."""
from chains import *
def run(L=4,bmax=6):
    one=0; pig=0; pig_classes=Counter()
    KB=U(0,0); legsB=[e for e in edges(KB) if cls(e)!=1]
    PB=[paths(e,KB,L,False) for e in legsB]
    for bA in [b for b in range(-bmax,bmax+1) if b!=0]:
        KA=U(0,bA); legsA=[e for e in edges(KA) if cls(e)!=1]
        PA=[paths(e,KA,L,False) for e in legsA]
        for p1 in PB[0]:
          t1=tris(p1)
          for p2 in PB[1]:
            t2=tris(p2)
            if t1&t2 or KA in t1|t2: continue
            for i in (0,1):
              for a1 in PA[i]:
                ta=tris(a1)
                if KB in ta: continue
                rh={frozenset(r[:2]):r for r in a1[0]}
                res=[]
                for pb,tb in ((p1,t1),(p2,t2)):
                    c=0; bad=None
                    for r in pb[0]:
                        ra=rh.get(frozenset(r[:2]))
                        if ra is None:
                            if r[0] in ta or r[1] in ta: bad=bad or 'tri'
                            continue
                        if cls(ra[2])==cls(r[2]): bad='sameclass'
                        else: c+=1
                    res.append((c,bad))
                if set(res)=={(0,None),(1,None)}: one+=1
                for j in (0,1):
                    if res[j]==(1,None) and res[1-j][1]=='sameclass':
                        pig+=1; pig_classes[cls(a1[1]) if not a1[0] else cls(a1[0][0][2])]+=1
    return one,pig,pig_classes
o,p,pc=run()
print('strict L=4: one-crossing pictures =',o,'; pigeonhole reached (1 good crossing + same-class collision with other half) =',p,'by A-class',dict(pc))
