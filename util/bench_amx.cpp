/**
 * @file bench_amx.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Benchmarks peak performance and latency of amx mmm instruction
 * 
 */

#include <cstdint>
#include <stdio.h>
#include <immintrin.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <papi.h>
#include "APIB/APIB_helper.h"

#define NUM_RUNS 100000
#define NUM_RUNS_WARMUP 10000
#define NUM_REPEATS 100

extern "C" uint64_t peak_func(uint64_t counter);
extern "C" uint64_t latency_func(uint64_t counter);

void handle_error (int retval) {
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

// takes function pointer to asm implementation of function as argument
double bench(uint64_t (*func)(uint64_t)) {
    // warmup
    uint64_t num_ins = func(NUM_RUNS_WARMUP);
    
    int retval; 
    int EventSet = PAPI_NULL;
    long_long values[1];

    // init library
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT)
        handle_error(retval);
    
    // event set
    retval = PAPI_create_eventset(&EventSet);
    if (retval != PAPI_OK)
        handle_error(retval);

    // add cycles to events
    retval = PAPI_add_event(EventSet, PAPI_TOT_CYC);
    if (retval != PAPI_OK)
        handle_error(retval);

    // start
    retval = PAPI_start(EventSet);
    if (retval != PAPI_OK)
        handle_error(retval);
    
    // run target code
    for (int i = 0; i < NUM_REPEATS; i++) {
        func(NUM_RUNS);
    }
    
    // stop
    retval = PAPI_stop(EventSet, values);
    if (retval != PAPI_OK)
        handle_error(retval);
    
    // return in unit instructions per cycle
    return ((double)num_ins * NUM_RUNS * NUM_REPEATS) / values[0];
}

int main() {
    // configure to use amx
    if (!set_tiledata_use()) exit(-1);
    __tilecfg tile_data = {0};
    // load tile config
    init_tile_config (&tile_data, MAX_ROWS, MAX_COLSB);

    // actual benchmarks
    printf("peak perf: %lf ops/cycle\n", bench(peak_func));
    // invert to get cycles / instruction for latency
    printf("latency: %lf cycles\n", 1.0 / bench(latency_func));

}
