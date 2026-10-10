from sweep2 import report
import sys
for (k,ell,n) in [(1,0,4),(1,0,5),(1,0,6),(1,1,6),(1,1,7),(1,2,8),(1,2,9),(1,3,10),(2,0,4),(2,0,5),(2,1,6)]:
    try:
        r=report(k,ell,n)
        print(f"k={r['k']} ell={r['ell']} w={r['w']} n={r['n']}: mono={r['nmono']:5d} | STRICT viol={len(r['strict']):5d} | WEAK viol={len(r['anysign']):5d}", flush=True)
        if r['strict']:  g,mo,d,A,B,T=r['strict'][0];  print(f"    worst STRICT: alpha={mo} coeff={d} A={A} B={B} TOT={T} deficit={g}", flush=True)
        if r['anysign']: g,mo,d,A,B,T=r['anysign'][0]; print(f"    worst WEAK  : alpha={mo} coeff={d} L1-TOT={g} TOT={T}", flush=True)
    except Exception as ex:
        print(f"k={k} ell={ell} n={n}: ERROR {type(ex).__name__} {ex}", flush=True)
