"""(2B*_d) and (2B**_d) as stated in Q150-c1 Convention 2.4, plus support analysis."""
from itertools import product

def words(n):
    return list(product((0,1), repeat=n))

def relations_2Bstar(d):
    """Yield (u1,v1,u2,v2) meaning sigma(u1)sigma(v1) = sigma(u2)sigma(v2)."""
    for delta in range(1, d):
        for pi in words(delta-1):
            for pip in words(d-delta-1):
                for sg in words(delta-1):
                    u0 = pi+(0,)+pip ; v0 = pip+(0,)+sg
                    u1 = pi+(1,)+pip ; v1 = pip+(1,)+sg
                    yield (u0,v0,u1,v1,delta)

def relations_2Bstarstar(d):
    for delta in range(1, d):
        for a in words(delta-1):
            for ap in words(d-delta-1):
                for b in words(d-delta-1):
                    for bp in words(delta-1):
                        yield (a+(0,)+ap, b+(0,)+bp, a+(1,)+ap, b+(1,)+bp, delta)

def support_ok(Z, rels):
    """Z a frozenset of words. Check [u0 in Z and v0 in Z] <=> [u1 in Z and v1 in Z]."""
    for u0,v0,u1,v1,_ in rels:
        if ((u0 in Z and v0 in Z) != (u1 in Z and v1 in Z)):
            return False
    return True
