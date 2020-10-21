/*
* NAME: Copyright (c) 2020, Biren Patel
* LISC: MIT License
* DESC: Unit and integration tests for random number generator library. Most of
* these tests are slow as monte carlo simulations are required in many cases to 
* approximately verify function behavior.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <immintrin.h>
#include <time.h>

#include "timeit.h"
#include "random.h"
#include "unity\unity.h"

/******************************************************************************/

void setUp(void) {}
void tearDown(void) {}

#define BIG_SIMULATION 2500000
#define MID_SIMULATION 500000
#define SMALL_SIMULATION 50000

/*******************************************************************************
Given a seed (the answer to life, the universe, everything), does PCG 64i output
the same stream on two different random_t variables?
*/

void test_deterministic_seed_pcg_output(void)
{
    //arrange
    random_t rng_1 = rng_init(42);
    assert(rng_1.state != 0 && "rdrand failure");
    
    random_t rng_2 = rng_init(42);
    assert(rng_2.state != 0 && "rdrand failure");
    
    //act-assert
    for (size_t i = 0; i < BIG_SIMULATION; i++)
    {
        uint64_t result_1 = rng_next(&rng_1);
        uint64_t result_2 = rng_next(&rng_2);
        
        TEST_ASSERT_EQUAL_UINT64(result_1, result_2);
    }
}

/*******************************************************************************
Check that rng_bias is correct by monte carlo simulation on probabilites of
1/256 through 255/256. At 1,000,000 simulations with floating precision, the
tolerance is set as +/- 0.0015.
*/

void test_monte_carlo_of_rng_bias_at_256_bits_of_resolution(void)
{
    //arrange
    struct tuple {float actual; float expected;} result[255];
    random_t rng = rng_init(0);
    assert(rng.state != 0 && "rdrand failure");
    float expected_counter = 0.0;
    
    //act
    for (size_t i = 0; i < 255; i++)
    {
        int success = 0;
        uint64_t numerator = i + 1;
        expected_counter += 0.00390625f;
        
        for (size_t j = 0; j < BIG_SIMULATION; j++)
        {
            if (rng_bias(&rng, numerator, 8) & 1) success++;
        }
        
        result[i].actual = ((float) success) / BIG_SIMULATION;
        result[i].expected = expected_counter;
    }
    
    //assert
    for (size_t i = 0; i < 255; ++i)
    {        
        TEST_ASSERT_FLOAT_WITHIN(.001f, result[i].expected, result[i].actual);
    }        
}

/*******************************************************************************
Given an input stream with bits biased to .125 probability of success, output
a stream of 135 bits with unbiased bits. The input stream has no autocorrelation
and 135 was selected as a non-multiple of 64.
*/

void test_von_neumann_debiaser_outputs_all_unbiased_bits(void)
{
    //arrange
    random_t rng = rng_init(0);
    assert(rng.state != 0 && "rdrand failure");
    
    uint64_t input_stream[35]; //2240 input bits should be enough
    uint64_t output_stream[3]; //space for 192 but only care about first 135
    float results[135] = {0}; //holds monte carlo result for each output bit
    
    stream_t info;
    
    //act-assert
    for (size_t i = 0; i < MID_SIMULATION; i++)
    {
        //populate the input stream with biased bits
        for (size_t j = 0; j < 35; ++j)
        {
            input_stream[j] = rng_bias(&rng, 32, 8);
        }
        
        //populate the output stream with unbiased bits
        info = rng_vndb(input_stream, output_stream, 2240, 135);
        TEST_ASSERT_EQUAL_INT(135, info.filled);
        
        //walk through the output stream and assess results of this round
        for (size_t k = 0; k < 135; k++)
        {
            if ((output_stream[k/64] >> (k % 64)) & 1) results[k]++;
        }
    }
    
    //assert
    for (size_t i = 0; i < 135; ++i)
    {
        results[i] /= MID_SIMULATION;
        TEST_ASSERT_FLOAT_WITHIN(.01f, 0.5f, results[i]);
    }
}

/*******************************************************************************
Given a time series of 1010101010...10 the autocorrelation at any lag k should
alternate between 1 and -1.
*/

void test_cyclic_autocorrelation_of_alternating_bitstream(void)
{
    //arrange
    double results[64] = {0};
    random_t rng = rng_init(0);
    assert(rng.state != 0 && "rdrand failure");
    
    uint64_t input_stream[100000] = {0};
    uint64_t word = 0xAAAAAAAAAAAAAAAA;
    
    for (size_t i = 0; i < 100000; i++)
    {
        input_stream[i] = word;
    }
    
    //act
    for (size_t i = 0; i < 64; i++)
    {
        results[i] = rng_cycc(input_stream, 6400000, i);
    }
    
    //assert
    for (size_t i = 0; i < 64; i++)
    {
        if (i & 1) TEST_ASSERT_EQUAL_FLOAT(-1, (float) results[i]);
        else TEST_ASSERT_EQUAL_FLOAT(1, (float) results[i]);
    }
}

/*******************************************************************************
Since the SIMD implemntation is quite tricky, I need to ensure that each 64 bit
block is actually genreated from an independent PCG stream over two steps. So,
this unit test generates those individual streams from single PCG32i streams,
and then checks those values against the corresponding vector block. Since PCG
32i is not a part of the API, I include its entire implementation here. Quick
and dirty stdlib rand for the four seeds for a more dynamic unit test.
*/

uint64_t mix (uint64_t value)
{
    value ^= value >> 30;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27;
    value *= 0x94d049bb133111ebULL;
    value ^= value >> 31;
    return value;
}

struct pcg32i
{
    uint32_t current;
    uint32_t increment;
};

uint32_t pcg32i_next(struct pcg32i *state)
{
    uint32_t x = state->current;
    state->current = state->current * 747796405U + state->increment;
    uint32_t fx = ((x >> ((x >> 28U) + 4U)) ^ x) * 277803737U;
    return (fx >> 22U) ^ fx;
}

void test_simd_pcg_32_bit_insecure_generator(void)
{
    //arrange
    srand((unsigned) time(NULL));
    uint64_t a = (uint64_t) rand();
    uint64_t b = (uint64_t) rand();
    uint64_t c = (uint64_t) rand();
    uint64_t d = (uint64_t) rand();
    
    simd_random_t simd_rng = simd_rng_init(a,b,c,d);
    
    struct pcg32i rng_1;
    struct pcg32i rng_2;
    struct pcg32i rng_3;
    struct pcg32i rng_4;
    
    rng_1.current = (uint32_t) mix(a);
    rng_1.increment = (uint32_t) (mix(mix(a)) | 1);
    
    rng_2.current = (uint32_t) mix(b);
    rng_2.increment = (uint32_t) (mix(mix(b)) | 1);
    
    rng_3.current = (uint32_t) mix(c);
    rng_3.increment = (uint32_t) (mix(mix(c)) | 1);
    
    rng_4.current = (uint32_t) mix(d);
    rng_4.increment = (uint32_t) (mix(mix(d)) | 1);
    
    __m256i simd_out_vec;
    uint32_t *simd_out;
    
    //act-assert
    for (size_t i = 0; i < BIG_SIMULATION; i++)
    {
        simd_out_vec = simd_rng_next(&simd_rng);
        simd_out = (uint32_t *) &simd_out_vec;
        
        TEST_ASSERT_EQUAL_UINT32(simd_out[0], pcg32i_next(&rng_1));
        TEST_ASSERT_EQUAL_UINT32(simd_out[1], pcg32i_next(&rng_1));
        
        TEST_ASSERT_EQUAL_UINT32(simd_out[2], pcg32i_next(&rng_2));
        TEST_ASSERT_EQUAL_UINT32(simd_out[3], pcg32i_next(&rng_2));
        
        TEST_ASSERT_EQUAL_UINT32(simd_out[4], pcg32i_next(&rng_3));
        TEST_ASSERT_EQUAL_UINT32(simd_out[5], pcg32i_next(&rng_3));
        
        TEST_ASSERT_EQUAL_UINT32(simd_out[6], pcg32i_next(&rng_4));
        TEST_ASSERT_EQUAL_UINT32(simd_out[7], pcg32i_next(&rng_4));
    }
}

/*******************************************************************************
Benchmarks on 1 million draws.
*/

#define loop for (size_t i = 0; i < 1000000; i++)

void speed_test(void)
{
    simd_random_t simd_rng = simd_rng_init(10, 20, 30, 40);
    random_t rng = rng_init(50);
    init_timeit();       
    puts("\n~~~~~ Speed Tests ~~~~~");
    
    //base PCG 64i generator
    start_timeit();
    loop { rng_next(&rng); }
    end_timeit();
    printf("PCG Generator (64 Bits): %llu us\n", result_timeit(MICROSECONDS));
    
    //base PCG 64i generator at 4 calls (256 bits)
    start_timeit();
    loop
    {
        rng_next(&rng); rng_next(&rng);
        rng_next(&rng); rng_next(&rng);
    }
    end_timeit();
    printf("PCG Generator (256 Bits): %llu us\n", result_timeit(MICROSECONDS));
    
    //SIMD generator at 1 call (256 bits)    
    start_timeit();
    loop 
    { simd_rng_next(&simd_rng); }
    end_timeit();
    printf("SIMD Generator (256 Bits): %llu us\n", result_timeit(MICROSECONDS));
    
    //rng bias at 8 generator calls
    start_timeit();
    loop { rng_bias(&rng, 1, 8); }
    end_timeit();
    printf("RNG Bias: %llu us\n", result_timeit(MICROSECONDS));
    
    //rng binomial at no additional generator calls (overhead only)
    start_timeit();
    loop { rng_bino(&rng, 64, 1, 8); }
    end_timeit();
    printf("RNG Binomial: %llu us\n", result_timeit(MICROSECONDS));
}

/******************************************************************************/

int main(void)
{
    UNITY_BEGIN();
        RUN_TEST(test_deterministic_seed_pcg_output);
        RUN_TEST(test_monte_carlo_of_rng_bias_at_256_bits_of_resolution);
        RUN_TEST(test_von_neumann_debiaser_outputs_all_unbiased_bits);
        RUN_TEST(test_cyclic_autocorrelation_of_alternating_bitstream);
        RUN_TEST(test_simd_pcg_32_bit_insecure_generator);
    UNITY_END();
    
    speed_test();
    
    return EXIT_SUCCESS;
}
