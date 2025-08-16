/**
 * @file eigen_bench_skeleton.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Skeleton to run benchmarks of eigen matrix multiplication
 * 
 */
 
#include <cstdint>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <math.h>
#include <Eigen/Dense>
#include "APIB/APIB_benchmark_specs.h"

using scalar_type = 
/*GEN SCALAR_TYPE*/ uint32_t
;
using matrix = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

// sizeof(data) * limit ≥ cache_size * associativity
// returns at least N_RUN_MIN
template<int n, int m, int p, int bits_A, int bits_B, int bits_C>
int get_limit() {
    int lim = (int)(ceil(((double)(CACHE_SIZE) * CACHE_ASSOC) / 
                (sizeof(scalar_type) * (n*p + p*m + n*m))));
    return lim >= N_RUN_MIN ? lim : N_RUN_MIN;
}

// init input
void fill_matrix(matrix& A, int rows, int cols) {
    for(int i=0; i < rows; i++) {
        for(int j=0; j < cols; j++) {
            A(i, j) = (scalar_type)(i+j);
        }
    }
}

matrix eigen_mmm(matrix& A, matrix& B, matrix& C, scalar_type alpha, scalar_type beta) {
    return alpha*C + beta*A*B;
}

void handle_error (int retval) {
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

// A: nxp, B: pxm, C: nxm
void benchmark(matrix* A, matrix* B, matrix* C, int n, int m, int p, int alpha, int beta, int limit) {
    // warmup cpu with last one
    auto res_warmup = eigen_mmm(A[limit - 1], B[limit - 1], C[limit - 1], alpha, beta);
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
            auto mmm_res = eigen_mmm(A[run], B[run], C[run], alpha, beta);
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
        /*GEN MAT_N*/ 1000
        ;
    // A: nxp, B: pxm, C: nxm
    const int m = 
        /*GEN MAT_M*/ 1000
        ;
    const int p = 
        /*GEN MAT_P*/ 1000
        ;
    // A: nxp, B: pxm, C: nxm
    const int bits_A = 
        /*GEN BITS_A*/ 3
        ;
    const int bits_B = 
        /*GEN BITS_B*/ 3
        ;
    const int bits_C = 
        /*GEN BITS_C*/ 3
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

    matrix A[limit];
    matrix B[limit];
    matrix C[limit];

    for (int i = 0; i < limit; i++) {
        A[i] = matrix(n ,p);
        fill_matrix(A[i], n, p);
        B[i] = matrix(p, m);
        fill_matrix(B[i], p, m);
        C[i] = matrix(n, m);
        fill_matrix(C[i], n, m);
    }

    benchmark(A, B, C, n, m, p, alpha, beta, limit);
    return 0;
}
