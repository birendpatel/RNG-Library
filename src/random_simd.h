/*
* NAME: Copyright (c) 2020, Biren Patel
* LISC: MIT License
* DESC: 256-Bit SIMD API (AVX/AVX2 Instruction Sets)
*/

#ifndef SIMD_RANDOM_H
#define SIMD_RANDOM_H

#include <stdint.h>
#include <immintrin.h>

/*******************************************************************************
* NAME: simd_random_t
* DESC: internal state of the default vectorized PRNG
* @ state : contains state of 4 streams in lower 32 bits of each 64 bit block
* @ increment : contain stream identifiers in lower 32 bits of each 64 bit block
*******************************************************************************/
typedef struct
{
    __m256i state;
    __m256i increment;
} simd_random_t;

/*******************************************************************************
* NAME: simd_rng_init
* DESC: initialize a variable of type simd_random_t
* OUTP: zero state and increment in return type indicates rdrand failure
* @ seed : If any seed is zero, then all four streams will be non-determinstic
*******************************************************************************/
simd_random_t simd_rng_init
(
    const uint64_t seed_1,
    const uint64_t seed_2,
    const uint64_t seed_3,
    const uint64_t seed_4
);

/*******************************************************************************
* NAME: simd_rng_next
* DESC: generate 256-Bit psuedo random numbers via the default PRNG
* OUTP: random numbers not guaranteed equal to the updated state parameters
*******************************************************************************/
__m256i simd_rng_next (simd_random_t * const rng);

#endif
