from gmpy2 import mpz, bit_mask
from tqdm import trange


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
    
    Mp = bit_mask(p)  # (1 << p) - 1
    s = mpz(4)
    
    for _ in range(p - 2):
        s = fast_mod(s*s - 2, Mp, p)
    
    return not s
