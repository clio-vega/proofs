import sys; sys.path.insert(0,'/home/clio/projects/scratch/20261021c2')
from g3 import *
print('=== FULL positive control: N=4, I={2,3}, m=2, type B, eps=1 -- ALL SIX root families')
print('    populated (|I|>=2, |I^c|>=2, some j notin I below I, some j notin I above I).')
N=4; I=(2,3); m=2; w=wvars(N)
fs=diag_factors(N,I,w,'B')
for lab,f in fs:
    print(f'    {lab:16s}  {str(f):30s}  ord_S = {ord_S(f,w,m,1)}')
print(f'    #factors = {len(fs)};  sum of ord_S = {sum(ord_S(f,w,m,1) for _,f in fs)}')
print(f'    ord_S Atilde^B_(I,I) from SYMMETRISATION = {ord_S(Atilde(N,w,I,I,X="B"),w,m,1)}   predicted 1')
print()
print('    and the same I with m=1 (m notin I), the ord=0 half:')
print(f'    ord_S Atilde^B_(I,I) from SYMMETRISATION = {ord_S(Atilde(N,w,I,I,X="B"),w,1,1)}   predicted 0')
