# RNG Library

This is a small library of psuedo-random number functions for non-cryptographic (and non-secure) applications such as statistics, machine learning, and simulations. It requires an x86-64 Intel/AMD processor which contain the x86 RDRAND intruction and the AVX2 instruction set.

A significant portion of the library is designed to perform random number manipulations on bitarray structures. Bitarrays can be processed in 64-bit (or 256-bit with AVX2) chunks by generating 64 (256) samples simultaneously for certain probability distributions. This library might be an appropriate replacement for your existing infrastructure you find yourself e.g. generating uniform doubles in (0, 1), converting them one at a time to a target distribution function, and then setting/clearing/testing each individual bit.

All API documentation is contained in the appropriate header files. SISD headers target standard 64-bit applications, SIMD will leverage AVX2. For Windows users, you cannot compile the code using MinGW64 due to a bug in the compiler for AVX2 intrinsics (current as of 2020 October 19). Clang is fine. 

# Demo

```C
#include "random.h"
#include <stdint.h>
#include <stddef.h>

int main(void)
{
    uint64_t target[100] = {0};

    //non-deterministic seeding with zero
    random_t rng = rng_init(0);
    
    //set each bit in the target to 1 with probability p = 32/2^8 = 0.125
    //in total, 300 calls to the underlying PRNG are made.
    //A U~(0,1) transformation would have required 6400 calls.
    for (size_t i = 0; i < 100; i++)
    {
        target |= rng.bias(&rng.state, 32, 8);
    }
    
    //for finer-grained p, the total calls approaches 6400 in the limit.
    //For p = 1/2^32 = 2.3E-10 then 3,200 calls are made.
    for (size_t i = 0; i < 100; i++)
    {
        target |= rng.bias(&rng.state, 1, 32);
    }
}
```

# Task List
- [ ] SIMD Bernoulli Sampling
- [ ] SIMD Unit Tests
- [ ] SISD Unit Tests
