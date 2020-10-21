/*
* NAME: Copyright (c) 2020, Biren Patel
* LISC: MIT License
* DESC: 64-Bit SISD API
*/

#ifndef SISD_RANDOM_H
#define SISD_RANDOM_H

#include <stdint.h>

/*******************************************************************************
* NAME: random_t
* DESC: internal state of the default PRNG
* @ current : state value used to generate PRNG values
* @ increment : stream identifier
*******************************************************************************/
typedef struct
{
    uint64_t state;
    uint64_t increment;
} random_t;
    
/*******************************************************************************
* NAME: stream_t
* DESC: return type for some of the functions which process bitstreams
* @ used : generally, the number of bits used from the input stream
* @ filled : generally, the number of bits written into the output stream
*******************************************************************************/
typedef struct
{
    uint64_t used;
    uint64_t filled;
} stream_t;

/*******************************************************************************
* NAME: rng_init
* DESC: initialize a variable of type random_t
* OUTP: type random_t where zero state and increment indicates rdrand failure.
* @ seed : set seed = 0 for non-deterministic seeding.
*******************************************************************************/
random_t rng_init(const uint64_t seed);

/*******************************************************************************
* NAME: rng_generator
* DESC: generate a psuedo random number via the default PRNG.
* OUTP: random number not guaranteed equal to the updated state parameter.
*******************************************************************************/
uint64_t rng_next(random_t * const rng);

/*******************************************************************************
* NAME: rng_rand
* DESC: generate an unbiased psuedo random number
* OUTP : 0 if null state. non-null state will be updated
* @ min : inclusive lower bound
* @ max : inclusive upper bound
*******************************************************************************/
uint64_t rng_rand(random_t * const rng, const uint64_t min, const uint64_t max);

/*******************************************************************************
* NAME: rng_bias
* DESC: simultaneous generation of 64 iid bernoulli trials 
* OUTP: 64-bit word where each bit has probability p = n/2^m of success
* NOTE: m limits the total calls to rng_generator, so smaller m is faster code
* @ n : nonzero numerator of probability, strictly less than 2^m
* @ m : nonzero base 2 exponent less than or equal to 64
*******************************************************************************/
uint64_t rng_bias(random_t * const rng, const uint64_t n, const int m);

/*******************************************************************************
* NAME: rng_vndb
* DESC: Von Neumann Debiaser for iid biased bits with zero autocorrelation
* OUTP: dest is filled with stream_t.filled bits, which used stream_t.used bits
* NOTE: dest is not guaranteed to be filled to capacity
* @ src : binary bit stream of length n bits
* @ dest : binary bit stream of length m bits
*******************************************************************************/
stream_t rng_vndb 
(
    const uint64_t * restrict src, 
    uint64_t * restrict dest, 
    const uint64_t n, 
    const uint64_t m
);

/*******************************************************************************
* NAME: rng_cycc
* DESC: calculate the cyclic autocorrelation of an n-bit binary bitstream
* OUTP: lag-k cyclic correlation in inclusive range [-1.0, 1.0]
* @ src : binary bit stream of length n bits
* @ k : autocorrelation lag not exceeding total bits n
*******************************************************************************/
double rng_cycc(const uint64_t *src, const uint64_t n, const uint64_t k);

/*******************************************************************************
* NAME: rng_binomial
* DESC: sample from a binomial distribution X~(k,p) where p = n/2^m
* OUTP: number of successful trials
* @ k : total trials
* @ n : nonzero numerator of probability, strictly less than 2^m
* @ m : nonzero base 2 exponent less than or equal to 64
*******************************************************************************/
uint64_t rng_bino(random_t * const rng, uint64_t k, const uint64_t n, const int m);

#endif
