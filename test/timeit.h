/*
* NAME: Copyright (c) 2020, Biren Patel
* DESC: Windows macros for speed benchmarking. I wouldn't say I really know the
* nuances of QPF just yet, but from what I've seen this seems good  enough for 
* some rough single-threaded tests.
*/

#ifndef TIMEIT_H
#define TIMEIT_H

#include <windows.h>

#define NANOSECONDS 1000000000
#define MICROSECONDS 1000000
#define MILLISECONDS 1000
#define SECONDS 1

#define init_timeit()                                                          \
        LARGE_INTEGER start;                                                   \
        LARGE_INTEGER end;                                                     \
        LARGE_INTEGER Hz;                                                      \
        QueryPerformanceFrequency(&Hz)                                         \
        
#define start_timeit()                                                         \
        QueryPerformanceCounter(&start)                                        \
        
#define end_timeit()                                                           \
        QueryPerformanceCounter(&end)                                          \
        
#define result_timeit(unit)                                                    \
        (((end.QuadPart - start.QuadPart) * unit) / Hz.QuadPart)               \

#endif
