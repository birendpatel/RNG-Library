/*
* NAME: Copyright (c) 2020, Biren Patel
* LISC: MIT License
* DESC: PRNG library utilities implementation
*/

#include "random_utils.h"

#include <immintrin.h>

/*******************************************************************************
Retry loop for RDRAND x86 instruction. Per the Intel documentation, we give up 
after ten failures.
*/

bool rdrand(uint64_t *x)
{
    for (size_t i = 0; i < 10; i++)
    {
        if (_rdrand64_step(x))
        {
            return true;
        }
    }
    
    return false;
}

/*******************************************************************************
This is used to mix a user-supplied seed, it is Sebastiano Vigna's version of
Java's SplittableRandom: http://xoshiro.di.unimi.it/splitmix64.c but since its
just a one-off mixing function I removed the state increment that Vigna employs.
http://prng.di.unimi.it/splitmix64.c
*/

uint64_t rng_hash 
(
    uint64_t value
)
{
    value ^=  value >> 30;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27;
    value *= 0x94d049bb133111ebULL;
    value ^= value >> 31;
    
    return value;
}
