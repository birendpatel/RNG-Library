/*
* NAME: Copyright (c) 2020, Biren Patel
* LISC: MIT License
* DESC: PRNG library implementation
*/

#include "random_sisd.h"
#include "random_utils.h"
#include "bitarray.h"

#include <string.h>
#include <assert.h>
#include <limits.h>

/*******************************************************************************
Since this is a non-crypto statistics library, I use rdrand instead of rdseed
because it is A) faster since it doesn't require a pass through an extrator for
full entropy and B) less or almost zero-prone to underflow. David Johnston
explains this well in 13.4 of "Random Number Generators". Per Intel docs, rdrand
is retried up to ten times per variable, hence the else clause goto fuckery. For
PCG, ensure the increment is odd.
*/

random_t rng_init
(
    const uint64_t seed
)
{
    random_t rng;
    
    if (seed != 0) 
    {
        rng.state = rng_hash(seed);
        rng.increment = rng_hash(rng_hash(seed));
    }
    else
    {
        if (rdrand(&rng.state) && rdrand(&rng.increment))
        {
            goto success;
        }

        rng.state = 0;
        rng.increment = 0;
        goto terminate;
    }
    
    success:
        rng.increment |= 1;
    
    terminate:
        return rng;
}

/*******************************************************************************
This function uses a virtual machine to interpret a portion of the bit pattern
in the numerator parameter as executable bitcode. I wrote a short essay at url
https://stackoverflow.com/questions/35795110/ (username Ollie) to demonstrate
the concepts using 256 bits of resolution.
*/

uint64_t rng_bias 
(
    random_t * const rng, 
    const uint64_t n, 
    const int m
)
{
    assert(rng != NULL && "generator is null");
    assert(n != 0 && "probability is 0");
    assert(m > 0 && m <= 64 && "invalid base 2 exponent");
    
    uint64_t accumulator = 0;
    
    for (int pc = __builtin_ctzll(n); pc < m; pc++)
    {
        switch ((n >> pc) & 1)
        {
            case 0:
                accumulator &= rng_next(rng);
                break;
                
            case 1:
                accumulator |= rng_next(rng);
                break;
        }
    }
    
    return accumulator;
}

/*******************************************************************************
Von Neumann Debiaser for biased bits with no autocorrelation. Feed a low entropy
n-bit bitstream into the debiaser, get a high-entropy at-most-m-bit bitstream.
It may be that not all source bits are used and/or not all destination bits are
filled. The source is read as consecutive bit-pairs. The destination must be
zeroed out before the main loop since bitwise-or is used to set the bit array.
*/

stream_t rng_vndb 
(
    const uint64_t * restrict src, 
    uint64_t * restrict dest, 
    const uint64_t n, 
    const uint64_t m
)
{  
    assert(src != NULL && "null source");
    assert(dest != NULL && "null dest");
    assert(n != 0 && "nothing to read");
    assert(m != 0 && "nowhere to write");
    assert(n % 2 == 0 && "cannot process odd-length bitstream");
    
    u64_bitarray(src);
    u64_bitarray(dest);
    
    uint64_t write_pos = 0;
    uint64_t read_pos = 0;
    
    stream_t info = {.used = 0, .filled = 0};
    memset(dest, 0, (m-1)/CHAR_BIT + 1);
        
    while (read_pos < n)
    {        
        switch (u64_bitarray_mask_at(src, read_pos, 0x3ULL))
        {
            case 1:
                u64_bitarray_set(dest, write_pos);
                write_pos++;
                break;
            case 2:
                write_pos++;
                break;
        }
        
        read_pos += 2;
        
        if (write_pos == m) goto dest_filled;
    }
    
    dest_filled:
        info.used = read_pos;
        info.filled = write_pos;
        return info;
}

/*******************************************************************************
Cyclic lag-K autocorrelation of an n-bit stream. This uses the SCC algorithm
from Donald Knuth as the base and adds the binary bit stream simplification
from David Johnston's "Random Number Generators". 
*/

double rng_cycc
(
    const uint64_t *src, 
    const uint64_t n, 
    const uint64_t k
)
{
    assert(src != NULL && "data pointer is null");
    assert(n != 0 && "no data");
    assert(k < n && "lag exceeds length of data");
    
    u64_bitarray(src);
    
    uint64_t i = 0;
    uint64_t x1 = 0;
    uint64_t x2 = 0;
    
    while (i < n)
    {
        if (u64_bitarray_test(src, i))
        {            
            if (u64_bitarray_test(src, (i + k) % n))
            {
                x1++;
            }
            
            x2++;
        }
        
        i++;
    }    
    
    double numerator = 
        ((double) n * (double) x1 - ((double) x2 * (double) x2));
    
    double denominator =
        ((double) n * (double) x2 - ((double) x2 * (double) x2));
    
    assert(numerator/denominator >= -1.0 && "lower bound violation");
    assert(numerator/denominator <= 1.0 && "upper bound violation");
    
    return numerator/denominator;
}

/*******************************************************************************
Bitmask rejection sampling technique that Apple uses in their 2008 arc4random C 
source. I made minor adjustments for a variable lower bound and inclusive upper 
bound. I also throw away the random number after failure instead of attempting 
to use the upper bits.
*/

uint64_t rng_rand
(
    random_t * const rng, 
    const uint64_t min, 
    const uint64_t max
)
{    
    assert(rng != NULL && "generator is null");
    assert(min < max && "bounds violation");
    
    uint64_t sample;
    uint64_t scaled_max = max - min;
    uint64_t bitmask = ~((uint64_t) 0) >> __builtin_clzll(scaled_max);
    
    assert(__builtin_clzll(bitmask) == __builtin_clzll(scaled_max) && "bad mask");
    assert(__builtin_popcountll(bitmask) == 64 - __builtin_clzll(scaled_max) && "bad mask");
    
    do 
    {
        sample = rng_next(rng) & bitmask;
    } 
    while (sample > scaled_max);
    
    assert(sample <= scaled_max && "scaled bounds violation");
    
    return sample + min;
}

/*******************************************************************************
Generate a number from a binomial distribution by simultaneous simulation of
64 iid bernoulli trials per loop.
*/

uint64_t rng_bino
(
    random_t * const rng, 
    uint64_t k, 
    const uint64_t n, 
    const int m
)
{
    assert(rng != NULL && "generator is null");
    assert(n != 0 && "probability is 0");
    assert(m > 0 && m <= 64 && "invalid base 2 exponent");
    assert(k != 0 && "no trials");
    
    uint64_t success = 0;
    
    for (; k > 64; k-= 64)
    {
        success += (uint64_t) __builtin_popcountll(rng_bias(rng, n, m));
    }
    
    return success + (uint64_t) __builtin_popcountll(rng_bias(rng, n, m) >> (64 - k));
}

/*******************************************************************************
The following code is originally Copyright 2014 Melissa O'Neill pcg_random.org,
Licensed under the Apache License, Version 2.0. 

This is the insecure 64-bit output PCG denoted as pcg_output_rxs_m_xs_64_64.
I have made a few changes. First, the seeding is not done with dev/urandom or
stack variable XORing which O'Neill uses. I have placed the seeding in rng_init
which uses either RDRAND or Vigna's modified SplitMix64 which is in random_utils
and located at http://prng.di.unimi.it/splitmix64.c . For non-deterministic
seeding, rdrand is used directly. For the actual generator, I have simply placed
all the relevant macros that O'Neill uses into a single function, with hex
constants and variable name changes for some simplicity. The LCG is not modified
and the permutation step is not modified.
*/

uint64_t rng_next
(
    random_t * const rng
)
{
    uint64_t x = rng->state;
    
    rng->state = rng->state * 0x5851F42D4C957F2DULL + rng->increment;
    
    uint64_t fx = ((x >> ((x >> 59ULL) + 5ULL)) ^ x) * 0xAEF17502108EF2D9ULL;
    
    return (fx >> 43ULL) ^ fx;
}
