/**
 * @file APIB_bench_skeleton.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Skeleton to run benchmarks of APIB matrix multiplication
 * 
 */
 
#include <cstdint>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <math.h>
#include "APIB/APIB.h"
#include "APIB/APIB_benchmark_specs.h"

// sizeof(data) * limit ≥ cache_size * associativity
// returns at least N_RUN_MIN
template<int n, int m, int p, int bits_A, int bits_B, int bits_C>
int get_limit() {
    int space_A = APIB_Mat<bits_A, n, p, uint64_t>::get_size(); 
    int space_B = APIB_Mat<bits_B, p, m, uint64_t>::get_size(); 
    int space_C = APIB_Mat<bits_C, n, m, uint64_t>::get_size();

    int lim = (int)(ceil(((double)(CACHE_SIZE) * CACHE_ASSOC) / 
                (space_A + space_B + space_C)));

    return lim >= N_RUN_MIN ? lim : N_RUN_MIN;
}

// init input
template<typename buf_type>
void fill_matrix(buf_type* A, int rows, int cols, int bits) {
    for(int i=0; i < rows; i++) {
        for(int j=0; j < cols; j++) {
            A[cols*i+j] = (buf_type)((i+j) % (1L << bits));
        }
    }
}

void handle_error (int retval) {
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

// A: nxp, B: pxm, C: nxm
template <int num_bits_alpha, int num_bits_beta, int num_bits_A, int num_bits_B, int num_bits_C, int n, int m, int p, typename type_A, typename type_B, typename type_C>
void benchmark(APIB_Mat<num_bits_A, n, p, type_A>* A, APIB_Mat<num_bits_B, p, m, type_B>* B, APIB_Mat<num_bits_C, n, m, type_C>* C, int alpha, int beta, int limit) {
    // warmup cpu with last one
    auto mmm_res_warmup = mmm<num_bits_alpha, num_bits_beta>(A[limit - 1], B[limit - 1], C[limit - 1], alpha, beta);
    int retval;
    int EventSet = PAPI_NULL;
    long_long values[2];

    // init papi
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT)
        handle_error(retval);
    
    // event set
    retval = PAPI_create_eventset(&EventSet);
    if (retval != PAPI_OK)
        handle_error(retval);

    // add total cycles
    retval = PAPI_add_event(EventSet, PAPI_TOT_CYC);
    if (retval != PAPI_OK)
        handle_error(retval);
    
    // add lvl 3 cache misses
    retval = PAPI_add_event(EventSet, PAPI_L3_TCM);
    if (retval != PAPI_OK)
        handle_error(retval);

    double work = limit*(2.0*n*m*p + 3.0*n*m);
    double perf = 0;
    double op_int = 0;
    for (int repeat = 0; repeat < N_REPEATS; repeat++) {
        // start measurement
        retval = PAPI_start(EventSet);
        if (retval != PAPI_OK)
            handle_error(retval);

        for (int run = 0; run < limit; run++) {
            auto mmm_res = mmm<num_bits_alpha, num_bits_beta>(A[run], B[run], C[run], alpha, beta);
        }

        // read from measurement
        retval = PAPI_read(EventSet, values);
        if (retval != PAPI_OK)
            handle_error(retval);

        // stop measurement
        retval = PAPI_stop(EventSet, values);
        if (retval != PAPI_OK)
            handle_error(retval);

        perf += work/values[0];
        op_int += work/(CACHE_BLOCK_SIZE * values[1]);

        // reset counters
        retval = PAPI_reset(EventSet);
        if (retval != PAPI_OK)
            handle_error(retval);
    }

    perf /= N_REPEATS;
    op_int /= N_REPEATS;

    /*GEN VERBOSE START*/ 
    printf("Performance: %lf ops/cycle\n", perf);
    printf("Operational Intensity: %lf ops/byte\n", op_int);
    /*GEN VERBOSE END*/
    FILE* res = fopen(
    /*GEN INSERT RESFILE*/ "res.txt"
        , "w");
    fprintf(res, "%lf\n", perf);
    fprintf(res, "%lf\n", op_int);
    fclose(res);
}

int main() {
    // A: nxp, B: pxm, C: nxm
    const int n = 
    /*GEN MAT_N*/ 2000
        ;
    const int m = 
    /*GEN MAT_M*/ 2000
    ;
    const int p = 
    /*GEN MAT_P*/ 2000
    ;
    // A: nxp, B: pxm, C: nxm
    const int bits_A = 
    /*GEN BITS_A*/ 2
    ;
    const int bits_B = 
    /*GEN BITS_B*/ 2
    ;
    const int bits_C = 
    /*GEN BITS_C*/ 2
    ;
    const int bits_alpha = 
    /*GEN BITS_ALPHA*/ 2
    ;
    const int bits_beta = 
    /*GEN BITS_BETA*/ 2
    ;

    srand(0);
    int alpha = rand() % (1L << bits_alpha);
    int beta = rand() % (1L << bits_beta);

    int limit = get_limit<n, m, p, bits_A, bits_B, bits_C>();

    uint64_t* A_data = (uint64_t*)malloc(n*p*sizeof(uint64_t));
    fill_matrix<uint64_t>(A_data, n, p, bits_A);
    uint64_t* B_data = (uint64_t*)malloc(p*m*sizeof(uint64_t));
    fill_matrix<uint64_t>(B_data, p, m, bits_B);
    uint64_t* C_data = (uint64_t*)malloc(n*m*sizeof(uint64_t));
    fill_matrix<uint64_t>(C_data, n, m, bits_C);
    
    APIB_Mat<bits_A, n, p, uint64_t> mat_A(A_data);
    APIB_Mat<bits_B, p, m, uint64_t> mat_B(B_data);
    APIB_Mat<bits_C, n, m, uint64_t> mat_C(C_data);

    free(A_data);
    free(B_data);
    free(C_data);
    
    APIB_Mat<bits_A, n, p, uint64_t> A[limit];
    APIB_Mat<bits_B, p, m, uint64_t> B[limit];
    APIB_Mat<bits_C, n, m, uint64_t> C[limit];

    for (int i = 0; i < limit; i++) {
        // creates deep copy
        A[i] = mat_A;
        B[i] = mat_B;
        C[i] = mat_C;
    }

    benchmark<bits_alpha, bits_beta>(A, B, C, alpha, beta, limit);
    return 0;
}
