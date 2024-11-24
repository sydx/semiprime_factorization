import fractions
import matplotlib.pyplot as plt
import numpy as np
import sympy

CURRENT_MAX_P = None
CACHED_PRIMES = None

def primes(max_p=10_000):
    if CURRENT_MAX_P is not None and CURRENT_MAX_P >= max_p:
        return [p for p in CACHED_PRIMES if p <= max_p]
    primes = []
    p = 1
    while p <= max_p:
        while not sympy.isprime(p): p += 1
        if p > max_p: continue
        primes.append(p)
        p += 1
    if CURRENT_MAX_P is None or max_p > CURRENT_MAX_P:
        CACHED_PRIMES = primes.copy()
    return primes

def top(max_p=10_000):
    prime_pairs = []
    p = 1
    for p in primes(max_p):
        if p > max_p: continue
        prime_pairs.append((p, p))
    return prime_pairs

def semi_top(order=1, max_p=10_000):
    if order == 0: return top(max_p=10_000)
    assert order >= 0
    prime_pairs = []
    for p in primes(max_p):
        prime_pairs_ = []
        for q in primes(p):
            if q == p: continue
            prime_pairs_.append((p, q))
        if len(prime_pairs_) > order - 1: prime_pairs.append(prime_pairs_[-order])
    return prime_pairs

def precompute_signatures(max_order=100, K=100, max_p=10_000):
    signatures = {}
    for order in range(0, max_order):
        print(f'Computing signature of order {order}')
        semi_top_of_order = semi_top(order, max_p=max_p)
        signatures[order] = {}
        for semiprime in [prime_pair[0] * prime_pair[1] for prime_pair in semi_top_of_order]:
            for modulo in range(2, max(K, order + 1)):
                if modulo not in signatures[order]:
                    signatures[order][modulo] = set()
                signatures[order][modulo].add(semiprime % modulo)
    return signatures

# This function implements a suboptimal search through the signatures.
# There are numerous way to optimize this search, but this is only a proof of concept.
def determine_order(semi_prime, signatures, min_order=-1):
    for order in sorted(signatures.keys()):
        if order < min_order: continue
        match = True
        for modulo in sorted(signatures[order].keys()):
            if (semi_prime % modulo) not in signatures[order][modulo]:
                match = False
                break
        if match: return order
    return None

def imbden(semi_prime):
    return int(np.ceil(np.sqrt(semi_prime)))

def factorize(semi_prime, signatures):
    min_order = -1
    while True:
        order = determine_order(semi_prime, signatures, min_order=min_order)
        if order is None: return None
        this_imbden = imbden(semi_prime)
        p, q = (this_imbden + order), (this_imbden - order)
        if p * q == semi_prime: return p, q
        min_order = order + 1
    return None
