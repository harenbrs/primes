from gmpy2 import mpz
from tqdm import trange


def mersenne(p):
    return (1 << mpz(p)) - 1


def fast_mod(k, Mp, p):
    if k < 0:
        k += Mp
    
    if k < Mp:
        return k
    
    if k == Mp:
        return 0
    
    return fast_mod((k & Mp) + (k >> p), Mp, p)


def lucas_lehmer(p, progress=True):
    if progress:
        range = trange
    
    Mp = mersenne(p)
    s = mpz(4)
    
    for _ in range(p - 2):
        s = fast_mod(s*s - 2, Mp, p)
    
    return not s
