/*
* NAME: Copyright (c) 2020, Biren Patel
* LISC: MIT License
* DESC: PRNG library for non-cryptographic non-secure purposes like statistics
* and simulations. This library depends on a modern Intel/AMD machine with the
* AVX2 instruction set. Include this header directly instead of the 64/256/utils
* components.
*/

#ifndef RANDOM_H
#define RANDOM_H

//General Utilities
#include "random_utils.h"

//64-Bit PRNG
#include "random_sisd.h"

//256-Bit PRNG
#include "random_simd.h"

#endif
