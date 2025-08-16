/**
 * @file mkl_bench_skeleton.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Skeleton to run benchmarks of mkl matrix multiplication
 * 
 */
 
#include <cstdint>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <math.h>
#include "mkl_cblas.h"
#include "mkl.h"
#include "APIB/APIB_benchmark_specs.h"

using op_type = 
/*GEN OP_TYPE*/ uint16_t
;

using res_type = uint32_t;

// sizeof(data) * limit ≥ cache_size * associativity
// returns at least N_RUN_MIN
template<int n, int m, int p, int bits_A, int bits_B, int bits_C>
int get_limit() {
    int lim = (int)(ceil(((double)(CACHE_SIZE) * CACHE_ASSOC) / 
                (sizeof(op_type) * (n*p + p*m) + sizeof(res_type) * n*m)));
    return lim >= N_RUN_MIN ? lim : N_RUN_MIN;
}


// init input
template <typename t>
void fill_matrix(t* A, int rows, int cols) {
    for(int i=0; i < rows; i++) {
        for(int j=0; j < cols; j++) {
            A[i*cols + j] = (t)(i+j);
        }
    }
}

void mkl_mmm(op_type* A, op_type* B, res_type* C, int n, int m, int p, int alpha, int beta) {
    int l = 0;
    /*GEN MKL_MMM*/ cblas_gemm_s16s16s32(CblasRowMajor, CblasNoTrans, CblasNoTrans, CblasFixOffset, n, m, p, (float)alpha, (short *)A, p, 0, (short *)B, m, 0, (float)beta, (int *)C, m, &l);
}

void handle_error (int retval) {
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

// A: nxp, B: pxm, C: nxm
void benchmark(op_type* A[], op_type* B[], res_type* C[], int n, int m, int p, int alpha, int beta, int limit) {
    // warmup cpu with last one
    mkl_mmm(A[limit - 1], B[limit - 1], C[limit - 1], n, m, p, alpha, beta);
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
            mkl_mmm(A[run], B[run], C[run], n, m, p, alpha, beta);
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
    const int n = 
    /*GEN MAT_N*/ 2000
    ;
    // A: nxp, B: pxm, C: nxm
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

    op_type* A[limit];
    op_type* B[limit];
    res_type* C[limit];

    for (int i = 0; i < limit; i++) {
        A[i] = (op_type*)mkl_malloc(n*p*sizeof(op_type), 64);
        fill_matrix<op_type>(A[i], n, p);
        B[i] = (op_type*)mkl_malloc(p*m*sizeof(op_type), 64);
        fill_matrix<op_type>(B[i], p, m);
        C[i] = (res_type*)mkl_malloc(n*m*sizeof(res_type), 64);
        fill_matrix<res_type>(C[i], n, m);
    }

    benchmark(A, B, C, n, m, p, alpha, beta, limit);

    return 0;
}
