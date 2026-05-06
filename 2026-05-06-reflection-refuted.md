# Reflection conjecture refuted (2026-05-06)

## The conjecture

At S_5, the σ_1-G1 atoms of the two non-trivial rank-2 shapes have multiplicity vectors (in the canonical basis B_4 = {q^c (1+q)^{4-2c}}) that are reversals of each other:

- A_{(3,2)} has B-vector (1, 5, 3), degree 4
- A_{(4,1)} has B-vector (1, 3, 5), degree 4

The conjecture: this reflection (3,2) ↔ (4,1), implemented by c → deg/2 − c on B-coordinates, is a structural involution on Specht modules that generalises beyond S_5.

## The n=6 data (refutation)

There are three rank-2, length-≤2 partitions of 6: (4,2), (3,3), (5,1). Their σ_1-G1 atom multiplicity vectors:

| shape   | A-degree | B-vector              |
|---------|----------|-----------------------|
| (4,2)   | 8        | (1, 5, 21, 30, 9)     |
| (3,3)   | 6        | (1, 6, 13, 6)         |
| (5,1)   | 6        | (1, 3, 9, 14)         |

For a "reflection pair" we need two shapes with the same A-degree whose B-vectors are reversals. Checking:

- (4,2) is alone at degree 8.
- (3,3) and (5,1) both have degree 6, but their B-vectors are (1,6,13,6) and (1,3,9,14), and neither is the reversal of the other (reverse of (1,6,13,6) is (6,13,6,1); of (1,3,9,14) is (14,9,3,1)).

So no two rank-2 length-≤2 shapes at S_6 have both the same A-degree *and* a vector-reversal relationship.

## Reproducer

`scratch/2026-05-06-twin-wgraph/n6_rank2_reflection.py` (output in `n6_rank2_output.txt`).

## Conclusion

The S_5 reflection (1,5,3) ↔ (1,3,5) is a 1-pair degree-coincidence, not a structural involution. There is no Specht-module involution behind it. Any further use of the "reflection" idea must be local to S_5 or based on different data.
