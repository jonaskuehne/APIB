/**
 * @file bandwidth_measure_bench_skeleton.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Skeleton to benchmark memory bandwidth of cpu
 * 
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <math.h>
#include <immintrin.h>
#include "mkl.h"
#include "APIB/APIB_benchmark_specs.h"

#define SCALAR_SIZE sizeof(double)

// sizeof(data) * limit ≥ cache_size * associativity
// returns at least 2
int get_limit(int n) {
  int lim = (int)(ceil((double)(CACHE_SIZE * CACHE_ASSOC) / (SCALAR_SIZE * (n*n + 2*n))));
  return lim >= N_RUN_MIN ? lim : N_RUN_MIN;
}

// init input
void fill_matrix(double A[], int rows, int cols) {
    for(int i=0; i < rows; i++) {
        for(int j=0; j < cols; j++) {
            A[cols*i+j] = (double) i / (j+1) + 1.0;
        }
    }
}

void mvm(double A[], double x[], double y[], int n) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 1.0, y, 1);
}

void handle_error (int retval) {
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

// A: nxp, B: pxm, C: nxm
void benchmark(double* A[], double* x[], double* y[], int n, int limit) {
    // warmup cpu with last one
    mvm(A[limit - 1], x[limit - 1], y[limit - 1], n);
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

    double bandwidth = 0;
    for (int repeat = 0; repeat < N_REPEATS; repeat++) {
        // start measurement
        retval = PAPI_start(EventSet);
        if (retval != PAPI_OK)
            handle_error(retval);

        for (int run = 0; run < limit; run++) {
            mvm(A[limit - 1], x[limit - 1], y[limit - 1], n);
        }

        // read from measurement
        retval = PAPI_read(EventSet, values);
        if (retval != PAPI_OK)
            handle_error(retval);

        // stop measurement
        retval = PAPI_stop(EventSet, values);
        if (retval != PAPI_OK)
            handle_error(retval);

        double memory_traffic = ((double)CACHE_BLOCK_SIZE)*values[1];
        bandwidth += memory_traffic / values[0];

        // reset counters
        retval = PAPI_reset(EventSet);
        if (retval != PAPI_OK)
            handle_error(retval);
    }

    bandwidth /= N_REPEATS;
    
    /*GEN VERBOSE START*/ 
    printf("bandwidth: %lf\n", bandwidth);
    /*GEN VERBOSE END*/
    FILE* res = fopen(
    /*GEN INSERT RESFILE*/ "res.txt"
        , "w");
    // write twice to be compatible with format in other benchmarks
    fprintf(res, "%lf\n", bandwidth);
    fprintf(res, "%lf\n", bandwidth);
    fclose(res);
}

int main() {
    // A: nxn, x: nx1, C: nx1
    const int n = 
    /*GEN MAT_N*/ 5000
        ;
    int limit = get_limit(n);
    double* A[limit];
    double* x[limit];
    double* y[limit];
    
    for (int i = 0; i < limit; i++) {
        A[i] = (double*)mkl_malloc(n*n*sizeof(double), CACHE_BLOCK_SIZE);
        fill_matrix(A[i], n, n);
        x[i] = (double*)mkl_malloc(n*sizeof(double), CACHE_BLOCK_SIZE);
        fill_matrix(x[i], n, 1);
        y[i] = (double*)mkl_malloc(n*sizeof(double), CACHE_BLOCK_SIZE);
        fill_matrix(y[i], n, 1);
    }
    
    benchmark(A, x, y, n, limit);

    for (int i = 0; i < limit; i++) {
        mkl_free(A[i]);
        mkl_free(x[i]);
        mkl_free(y[i]);
    }
    return 0;
}
