"""T0. The engine agrees with the head-sum definition. Positive control for everything below."""
import random
from fractions import Fraction as F
from eng import psi_pair, psi_head, admissible

random.seed(1024_2)
ok = bad = 0; fails = []
for n in range(1, 6):
    for r in range(1, min(n, 4)+1):
        for _ in range(25):
            while True:
                w = [F(random.randint(-9, 9), random.randint(1, 4)) for _ in range(n)]
                if admissible(w): break
            theta = [random.randint(-2, 3) for _ in range(r)]
            try:
                a = psi_pair(w, theta); b = psi_head(w, theta)
            except ZeroDivisionError:
                continue
            if a == b: ok += 1
            else:
                bad += 1
                if len(fails) < 3: fails.append((n, r, w, theta, a, b))
print(f"T0 psi_pair == psi_head : {ok}/{ok+bad}   (n<=5, r<=4, signs allowed, no hygiene)")
for f in fails: print("  FAIL", f)
