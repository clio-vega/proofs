"""control for wakeA.window: with RA = RJ = R* (class-0 b rhombus) it must reproduce shpast's 21 feasible tiles."""
import sys
sys.argv = ['wakeA.py', '6']
src = open('wakeA.py').read().split('\nnr = lambda')[0]
exec(src)
Rstar = CF['B'][(k0, 'b')][0]
for T, Tp in (('B', 'C'), ('C', 'B')): print(T, window(T, Tp, Rstar, Rstar))
