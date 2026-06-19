"""Lemma N4 (c=4 interior Number Lemma) + the two j in {8,10} finite residue checks.
Proves the §6 residual of (a,b,4) interior. 2026-06-19."""
import numpy as np

def Hval(a,b,j):
    return (a**3*b**3+9*a**3*b**2+26*a**3*b+24*a**3+12*a**2*b**3-4*a**2*b**2*j**2
     -8*a**2*b**2*j+108*a**2*b**2-12*a**2*b*j**2-48*a**2*b*j+312*a**2*b-8*a**2*j**2
     -64*a**2*j+288*a**2+47*a*b**3-20*a*b**2*j**2-64*a*b**2*j+423*a*b**2+6*a*b*j**4
     -12*a*b*j**3-30*a*b*j**2-384*a*b*j+1222*a*b+24*a*j**3-40*a*j**2-488*a*j+1128*a
     +60*b**3-24*b**2*j**2-120*b**2*j+540*b**2+6*b*j**4+12*b*j**3-42*b*j**2-696*b*j
     +1560*b-4*j**6+48*j**5-244*j**4+648*j**3-664*j**2-648*j+1440)

# ---- Lemma N4: 16 | H for all a==b (mod 2), all j (complete residue check mod 16) ----
bad=sum(1 for a in range(16) for b in range(a%2,16,2) for j in range(16) if Hval(a,b,j)%16)
print(f"Lemma N4:  16 | H  over residues mod 16, a==b mod2 -> failures = {bad}  "
      f"({'PROVEN' if bad==0 else 'FALSE'})")
sharp = any(Hval(a,b,j)%32 for a in range(32) for b in range(a%2,32,2) for j in range(32))
print(f"           sharp (32 fails somewhere): {sharp};  v2H(8,8,8) = "
      f"{(lambda n: (n&-n).bit_length()-1)(Hval(8,8,8))}")

# ---- j in {8,10}:  2^k_j | (a+2)_{j-3}(b+1)_{j-3} H(j)  (complete residue check mod 2^k) ----
def Hmod(a,b,j,M):
    a%=M;b%=M;j%=M
    def mul(*xs):
        r=np.ones_like(b)
        for x in xs:r=(r*(np.int64(x)%M))%M
        return r
    terms=[mul(a,a,a,b,b,b),9*mul(a,a,a,b,b),26*mul(a,a,a,b),24*mul(a,a,a),12*mul(a,a,b,b,b),
     -4*mul(a,a,b,b)*(j*j),-8*mul(a,a,b,b)*j,108*mul(a,a,b,b),-12*mul(a,a,b)*(j*j),-48*mul(a,a,b)*j,
     312*mul(a,a,b),-8*mul(a,a)*(j*j),-64*mul(a,a)*j,288*mul(a,a),47*mul(a,b,b,b),-20*mul(a,b,b)*(j*j),
     -64*mul(a,b,b)*j,423*mul(a,b,b),6*mul(a,b)*(j**4),-12*mul(a,b)*(j**3),-30*mul(a,b)*(j*j),
     -384*mul(a,b)*j,1222*mul(a,b),24*a*(j**3),-40*a*(j*j),-488*a*j,1128*a,60*mul(b,b,b),
     -24*mul(b,b)*(j*j),-120*mul(b,b)*j,540*mul(b,b),6*b*(j**4),12*b*(j**3),-42*b*(j*j),-696*b*j,
     1560*b,-4*(j**6),48*(j**5),-244*(j**4),648*(j**3),-664*(j*j),-648*j,1440]
    s=np.zeros_like(b)
    for t in terms:s=(s+(np.int64(t)%M))%M
    return s

def residue_check(j,k):
    M=1<<k;fails=0
    for a in range(M):
        b=np.arange(a%2,M,2,dtype=np.int64)
        fa=1
        for i in range(j-3):fa=(fa*((a+2-i)%M))%M
        fb=np.ones_like(b)
        for i in range(j-3):fb=(fb*((b+1-i)%M))%M
        prod=(((fa*fb)%M)*Hmod(a,b,j,M))%M
        fails+=int(np.count_nonzero(prod))
    print(f"j={j}: 2^{k} | (a+2)_{{{j-3}}}(b+1)_{{{j-3}}}H({j}) over a,b mod {M}, a==b mod2 "
          f"-> failures={fails}  ({'PROVEN' if fails==0 else 'FALSE'})")

residue_check(8,12)
residue_check(10,14)   # ~2 min
