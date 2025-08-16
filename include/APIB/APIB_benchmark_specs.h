/**
 * @file APIB_benchmark_specs.h
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains parameters for benchmarks
 * 
 */

#pragma once

#define CACHE_BLOCK_SIZE 64
#define CACHE_ASSOC 16
#define CACHE_SIZE (30*(1<<20))
#define N_RUN_MIN 5
#define N_REPEATS 5
