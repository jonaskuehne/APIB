/**
 * @file LUT_bench_skeleton.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Skeleton to run benchmarks of matrix multiplication using LUTs
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
#include "APIB/LUT.h"

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
template <typename align_type, int N_U, int M_U, int P_U, int n_u, int m_u, int p_u, int num_bits_A, int num_bits_B, int num_bits_C, int n, int m, int p, typename type_A, typename type_B, typename type_C>
void benchmark(APIB_Mat<num_bits_A, n, p, type_A>* A, APIB_Mat<num_bits_B, p, m, type_B>* B, APIB_Mat<num_bits_C, n, m, type_C>* C, int limit, align_type* lut) {
    // warmup cpu with last one
    auto mmm_res_warmup = lut_mmm<align_type>(A[limit - 1], B[limit - 1], C[limit - 1], N_U, M_U, P_U, n_u, m_u, p_u, lut);
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
            auto mmm_res = lut_mmm<align_type>(A[run], B[run], C[run], N_U, M_U, P_U, n_u, m_u, p_u, lut);
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
    /*GEN MAT_N*/ 512
        ;
    const int m = 
    /*GEN MAT_M*/ 512
    ;
    const int p = 
    /*GEN MAT_P*/ 512
    ;
    // A: nxp, B: pxm, C: nxm
    const int num_bits = 1;
    const int bits_A = num_bits;
    const int bits_B = num_bits;
    const int bits_C = num_bits;

    // LUT tile size
    const int n_u = 4;
    const int m_u = 4;
    const int p_u = 2;

    // block sizes
    // tried different sizes, does not seem to make a difference
    const int N_U = 128;
    const int M_U = 128;
    const int P_U = 128;

    srand(0);

    int limit = get_limit<n, m, p, bits_A, bits_B, bits_C>();

    using align_type = uint32_t;
    using mat_type = uint8_t;

    mat_type* A_data = (mat_type*)malloc(n*p*sizeof(mat_type));
    fill_matrix<mat_type>(A_data, n, p, bits_A);
    mat_type* B_data = (mat_type*)malloc(p*m*sizeof(mat_type));
    fill_matrix<mat_type>(B_data, p, m, bits_B);
    mat_type* C_data = (mat_type*)malloc(n*m*sizeof(mat_type));
    fill_matrix<mat_type>(C_data, n, m, bits_C);
    
    APIB_Mat<bits_A, n, p, mat_type> mat_A(A_data);
    APIB_Mat<bits_B, p, m, mat_type> mat_B(B_data);
    APIB_Mat<bits_C, n, m, mat_type> mat_C(C_data);

    free(A_data);
    free(B_data);
    free(C_data);
    
    APIB_Mat<bits_A, n, p, mat_type> A[limit];
    APIB_Mat<bits_B, p, m, mat_type> B[limit];
    APIB_Mat<bits_C, n, m, mat_type> C[limit];

    align_type* lut = generate_lut<align_type>(n_u, m_u, p_u, num_bits);

    for (int i = 0; i < limit; i++) {
        A[i] = mat_A;
        B[i] = mat_B;
        C[i] = mat_C;
    }

    benchmark<align_type, N_U, M_U, P_U, n_u, m_u, p_u>(A, B, C, limit, lut);
    return 0;
}
