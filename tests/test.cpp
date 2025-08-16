/**
 * @file test.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Test suite for APIB project
 * 
 */

#include "APIB/APIB.h"
#include "APIB/APIB_unpacking.h"
#include "APIB/APIB_scaled_unpacking.h"
#include "APIB/LUT.h"
#include <cstdint>
#include <cstdlib>

// returns the smallest type this number of bits fit into
template<int num_bits>
struct get_test_buf_type {
    using type = typename std::conditional<num_bits <= 8, uint8_t, 
          typename std::conditional<num_bits <= 16, uint16_t, 
          typename std::conditional<num_bits <= 32, uint32_t, 
          uint64_t>::type>::type>::type;
};

// fills array with pseudorandom numbers of specified number of bits
template<typename t>
void fill(t* buf, int rows, int cols, int bits) {
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            if (bits >= 64) {
                buf[row*cols + col] = (t)rand();
            } else if (bits >= 63) {
	            buf[row*cols + col] = (t)(rand() & (~(1L << bits)));
	    } else {
                buf[row*cols + col] = (t)(rand() & ((1L << bits) - 1));
            }
        }
    }
}

// packing tests per type
template<typename test_type>
void do_packing() {
    const int n = 1000;
    const int bits_type = 8 * sizeof(test_type);
    test_type* A = (test_type*)malloc(n*n*sizeof(test_type));
    test_type* check = (test_type*)malloc(n*n*sizeof(test_type));

    // first check with full type for edge cases
    fill<test_type>(A, n, n, bits_type);
    APIB_Mat<bits_type, n, n, test_type> M1(A);
    M1.template unpack_data<test_type>(check, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((test_type)check[i*n + j] != (test_type)A[i*n + j]) {
                printf("failed\n");
                exit(1);
            }
        }
    }

    // non-full type
    fill<test_type>(A, n, n, (bits_type >> 1) - 1);
    APIB_Mat<(bits_type >> 1) - 1, n, n, test_type> M2(A);
    M2.template unpack_data<test_type>(check, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((test_type)check[i*n + j] != (test_type)A[i*n + j]) {
                printf("failed\n");
                exit(1);
            }
        }
    }

    printf("passed\n");

    free(A);
    free(check);
}

// test packing for all relevant buffer types
void test_packing() {
    printf("running packing tests\n");
    printf("uint_8: ");
    do_packing<uint8_t>();
    printf("uint_16: ");
    do_packing<uint16_t>();
    printf("uint_32: ");
    do_packing<uint32_t>();
    printf("uint_64: ");
    do_packing<uint64_t>();
    printf("done with packing tests\n\n");
}

// unpacking tests per type
template<typename test_type>
void do_unpacking() {
    const int n = 1000;
    const int bits_type = 8 * sizeof(test_type);
    test_type* A = (test_type*)malloc(n*n*sizeof(test_type));

    // use full type to test for edge cases
    fill<test_type>(A, n, n, bits_type);
    APIB_Mat<bits_type, n, n, test_type> M1(A);

    test_type* check1 = (test_type*)malloc(n*M1.block_width*sizeof(test_type));

    unpack<test_type, typename get_base_type<bits_type>::type>(check1, M1.packed_data, M1.get_loc(0, 0), n, M1.block_width, M1.block_width, bits_type);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < M1.block_width; j++) {
            if ((test_type)check1[i*M1.block_width + j] != (test_type)A[i*n + j]) {
                printf("failed\n");
                exit(1);
            } 
        }
    }

    // non-full type
    fill<test_type>(A, n, n, (bits_type >> 1) - 1);
    APIB_Mat<(bits_type >> 1) - 1, n, n, test_type> M2(A);
    test_type* check2 = (test_type*)malloc((n-2)*M2.block_width*sizeof(test_type));

    unpack<test_type, typename get_base_type<(bits_type >> 1) - 1>::type>(check2, M2.packed_data, M2.get_loc(1, 1), n-2, M2.block_width-2, M2.block_width, (bits_type >> 1) - 1);

    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < M2.block_width - 2; j++) {
            if ((test_type)check2[i*M2.block_width + j] != (test_type)A[(i+1)*n + (j+1)]) {
                printf("failed\n");
                exit(1);
            } 
        }
    }
    printf("passed\n");

    free(A);
    free(check1);
    free(check2);
}

// unpacking tests for all relevant types
void test_unpacking() {
    printf("running unpacking tests\n");
    printf("uint_8: ");
    do_unpacking<uint8_t>();
    printf("uint_16: ");
    do_unpacking<uint16_t>();
    printf("uint_32: ");
    do_unpacking<uint32_t>();
    printf("uint_64: ");
    do_unpacking<uint64_t>();
    printf("done with unpacking tests\n\n");
}

// unpacking tests with scaling per type
template<typename test_type>
void do_scale_unpacking() {
    const int n = 1000;
    const int bits_type = 8 * sizeof(test_type);
    test_type* A = (test_type*)malloc(n*n*sizeof(test_type));

    // use full type for edge cases
    fill<test_type>(A, n, n, bits_type - 1);
    APIB_Mat<bits_type - 1, n, n, test_type> M1(A);

    test_type* check1 = (test_type*)malloc(n*M1.block_width*sizeof(test_type));

    scaled_unpack<test_type, typename get_base_type<bits_type - 1>::type>(check1, M1.packed_data, 2, M1.get_loc(0, 0), n, M1.block_width, M1.block_width, bits_type - 1);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < M1.block_width; j++) {
            if ((test_type)check1[i*M1.block_width + j] != (test_type)(2 * A[i*n + j])) {
                printf("failed\n");
                exit(1);
            } 
        }
    }

    // use non-full types
    fill<test_type>(A, n, n, (bits_type >> 1) - 1);
    APIB_Mat<(bits_type >> 1) - 1, n, n, test_type> M2(A);
    test_type* check2 = (test_type*)malloc((n-2)*M2.block_width*sizeof(test_type));

    scaled_unpack<test_type, typename get_base_type<(bits_type >> 1) - 1>::type>(check2, M2.packed_data, 3, M2.get_loc(1, 1), n-2, M2.block_width-2, M2.block_width, (bits_type >> 1) - 1);

    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < M2.block_width - 2; j++) {
            if ((test_type)check2[i*M2.block_width + j] != (test_type)(3 * A[(i+1)*n + (j+1)])) {
                printf("failed\n");
                exit(1);
            } 
        }
    }
    printf("passed\n");

    free(A);
    free(check1);
    free(check2);
}

// run tests for unpacking with scaling for all relevant types
void test_scale_unpacking() {
    printf("running scale unpacking tests\n");
    printf("uint_8: ");
    do_scale_unpacking<uint8_t>();
    printf("uint_16: ");
    do_scale_unpacking<uint16_t>();
    printf("uint_32: ");
    do_scale_unpacking<uint32_t>();
    printf("uint_64: ");
    do_scale_unpacking<uint64_t>();
    printf("done with scale unpacking tests\n\n");
}

// execute generic mmm test per config
template<int n, int m, int p, int alpha, int beta, int num_bits_A, int num_bits_B, int num_bits_C>
void do_mmm() {
    // use 64bit here to make sure the specified number of bits fits
    uint64_t* A = (uint64_t*)malloc(n*p*sizeof(uint64_t));
    uint64_t* B = (uint64_t*)malloc(p*m*sizeof(uint64_t));
    uint64_t* C = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    uint64_t* D_unpacked = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    fill<uint64_t>(A, n, p, num_bits_A);
    fill<uint64_t>(B, p, m, num_bits_B);
    fill<uint64_t>(C, n, m, num_bits_C);

    APIB_Mat<num_bits_A, n, p, uint64_t> mat_A(A);
    APIB_Mat<num_bits_B, p, m, uint64_t> mat_B(B);
    APIB_Mat<num_bits_C, n, m, uint64_t> mat_C(C);

    auto D = mmm<Ceil_Log2<alpha>::value, Ceil_Log2<beta>::value>(mat_A, mat_B, mat_C, alpha, beta);

    D.template unpack_data<uint64_t>(D_unpacked, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            uint64_t ref = (uint64_t)(alpha*C[i*m + j]);

            for (int k = 0; k < p; k++) {
                ref += (uint64_t)(beta*A[i*p + k] * B[k*m + j]);
            }

            if (D_unpacked[i*m + j] != ref) {
                printf("failed\n");
                exit(1);
            }
        }
    }

    printf("passed\n");
    free(A);
    free(B);
    free(C);
    free(D_unpacked);
}

// execute mmm test where the A matrix has two slices
template<int n_1, int n_2, int m, int p, int alpha, int beta, int num_bits_A_1, int num_bits_A_2, int num_bits_B, int num_bits_C>
void do_slice_mmm_2() {
    const int n = n_1 + n_2;

    // use smallest type to enable testing when the buffer types of the slices are different
    using type_a_1 = typename get_test_buf_type<num_bits_A_1>::type;
    using type_a_2 = typename get_test_buf_type<num_bits_A_2>::type;

    type_a_1* A_1 = (type_a_1*)malloc(n_1*p*sizeof(type_a_1));
    type_a_2* A_2 = (type_a_2*)malloc(n_2*p*sizeof(type_a_2));
    uint64_t* B = (uint64_t*)malloc(p*m*sizeof(uint64_t));
    uint64_t* C = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    uint64_t* D_unpacked = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    fill<type_a_1>(A_1, n_1, p, num_bits_A_1);
    fill<type_a_2>(A_2, n_2, p, num_bits_A_2);
    fill<uint64_t>(B, p, m, num_bits_B);
    fill<uint64_t>(C, n, m, num_bits_C);

    using S1 = APIB_Slice_Spec<num_bits_A_1, n_1, type_a_1>;
    using S2 = APIB_Slice_Spec<num_bits_A_2, n_2, type_a_2>;

    APIB_Slice_Mat<p, S1, S2> mat_A(A_1, A_2);
    APIB_Mat<num_bits_B, p, m, uint64_t> mat_B(B);
    APIB_Mat<num_bits_C, n, m, uint64_t> mat_C(C);

    auto D = mmm<Ceil_Log2<alpha>::value, Ceil_Log2<beta>::value>(mat_A, mat_B, mat_C, alpha, beta);

    D.template unpack_data<uint64_t>(D_unpacked, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            uint64_t ref = (uint64_t)(alpha)*C[i*m + j];

            for (int k = 0; k < p; k++) {
                if (i < n_1) {
                    ref += (uint64_t)(beta)*A_1[i*p + k] * B[k*m + j];
                } else {
                    ref += (uint64_t)(beta)*A_2[(i-n_1)*p + k] * B[k*m + j];
                }
            }

            if (D_unpacked[i*m + j] != ref) {
                printf("failed\n");
                exit(1);
            }
        }

    }

    printf("passed\n");

    free(A_1);
    free(A_2);
    free(B);
    free(C);
    free(D_unpacked);
}

// execute mmm test where the A matrix has three slices
template<int n_1, int n_2, int n_3, int m, int p, int alpha, int beta, int num_bits_A_1, int num_bits_A_2, int num_bits_A_3, int num_bits_B, int num_bits_C>
void do_slice_mmm_3() {
    const int n = n_1 + n_2 + n_3;

    // use smallest type to enable testing when the buffer types of the slices are different
    using type_a_1 = typename get_test_buf_type<num_bits_A_1>::type;
    using type_a_2 = typename get_test_buf_type<num_bits_A_2>::type;
    using type_a_3 = typename get_test_buf_type<num_bits_A_3>::type;

    type_a_1* A_1 = (type_a_1*)malloc(n_1*p*sizeof(type_a_1));
    type_a_2* A_2 = (type_a_2*)malloc(n_2*p*sizeof(type_a_2));
    type_a_3* A_3 = (type_a_3*)malloc(n_3*p*sizeof(type_a_3));
    uint64_t* B = (uint64_t*)malloc(p*m*sizeof(uint64_t));
    uint64_t* C = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    uint64_t* D_unpacked = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    fill<type_a_1>(A_1, n_1, p, num_bits_A_1);
    fill<type_a_2>(A_2, n_2, p, num_bits_A_2);
    fill<type_a_3>(A_3, n_3, p, num_bits_A_3);
    fill<uint64_t>(B, p, m, num_bits_B);
    fill<uint64_t>(C, n, m, num_bits_C);

    using S1 = APIB_Slice_Spec<num_bits_A_1, n_1, type_a_1>;
    using S2 = APIB_Slice_Spec<num_bits_A_2, n_2, type_a_2>;
    using S3 = APIB_Slice_Spec<num_bits_A_3, n_3, type_a_3>;

    APIB_Slice_Mat<p, S1, S2, S3> mat_A(A_1, A_2, A_3);
    APIB_Mat<num_bits_B, p, m, uint64_t> mat_B(B);
    APIB_Mat<num_bits_C, n, m, uint64_t> mat_C(C);

    auto D = mmm<Ceil_Log2<alpha>::value, Ceil_Log2<beta>::value>(mat_A, mat_B, mat_C, alpha, beta);

    D.template unpack_data<uint64_t>(D_unpacked, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            uint64_t ref = (uint64_t)(alpha)*C[i*m + j];

            for (int k = 0; k < p; k++) {
                if (i < n_1) {
                    ref += (uint64_t)(beta)*A_1[i*p + k] * B[k*m + j];
                } else if (i < n_1 + n_2) {
                    ref += (uint64_t)(beta)*A_2[(i-n_1)*p + k] * B[k*m + j];
                } else {
                    ref += (uint64_t)(beta)*A_3[(i-n_1-n_2)*p + k] * B[k*m + j];
                }
            }

            if (D_unpacked[i*m + j] != ref) {
                printf("failed\n");
                exit(1);
            }
        }

    }

    printf("passed\n");

    free(A_1);
    free(A_2);
    free(A_3);
    free(B);
    free(C);
    free(D_unpacked);
}

// execute mmm test where the A matrix has four slices
template<int n_1, int n_2, int n_3, int n_4, int m, int p, int alpha, int beta, int num_bits_A_1, int num_bits_A_2, int num_bits_A_3, int num_bits_A_4, int num_bits_B, int num_bits_C>
void do_slice_mmm_4() {
    const int n = n_1 + n_2 + n_3 + n_4;

    // use smallest type to enable testing when the buffer types of the slices are different
    using type_a_1 = typename get_test_buf_type<num_bits_A_1>::type;
    using type_a_2 = typename get_test_buf_type<num_bits_A_2>::type;
    using type_a_3 = typename get_test_buf_type<num_bits_A_3>::type;
    using type_a_4 = typename get_test_buf_type<num_bits_A_4>::type;

    type_a_1* A_1 = (type_a_1*)malloc(n_1*p*sizeof(type_a_1));
    type_a_2* A_2 = (type_a_2*)malloc(n_2*p*sizeof(type_a_2));
    type_a_3* A_3 = (type_a_3*)malloc(n_3*p*sizeof(type_a_3));
    type_a_4* A_4 = (type_a_4*)malloc(n_3*p*sizeof(type_a_4));
    uint64_t* B = (uint64_t*)malloc(p*m*sizeof(uint64_t));
    uint64_t* C = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    uint64_t* D_unpacked = (uint64_t*)malloc(n*m*sizeof(uint64_t));

    fill<type_a_1>(A_1, n_1, p, num_bits_A_1);
    fill<type_a_2>(A_2, n_2, p, num_bits_A_2);
    fill<type_a_3>(A_3, n_3, p, num_bits_A_3);
    fill<type_a_4>(A_4, n_4, p, num_bits_A_4);
    fill<uint64_t>(B, p, m, num_bits_B);
    fill<uint64_t>(C, n, m, num_bits_C);

    using S1 = APIB_Slice_Spec<num_bits_A_1, n_1, type_a_1>;
    using S2 = APIB_Slice_Spec<num_bits_A_2, n_2, type_a_2>;
    using S3 = APIB_Slice_Spec<num_bits_A_3, n_3, type_a_3>;
    using S4 = APIB_Slice_Spec<num_bits_A_4, n_4, type_a_4>;

    APIB_Slice_Mat<p, S1, S2, S3, S4> mat_A(A_1, A_2, A_3, A_4);
    APIB_Mat<num_bits_B, p, m, uint64_t> mat_B(B);
    APIB_Mat<num_bits_C, n, m, uint64_t> mat_C(C);

    auto D = mmm<Ceil_Log2<alpha>::value, Ceil_Log2<beta>::value>(mat_A, mat_B, mat_C, alpha, beta);

    D.template unpack_data<uint64_t>(D_unpacked, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            uint64_t ref = (uint64_t)(alpha)*C[i*m + j];

            for (int k = 0; k < p; k++) {
                if (i < n_1) {
                    ref += (uint64_t)(beta)*A_1[i*p + k] * B[k*m + j];
                } else if (i < n_1 + n_2) {
                    ref += (uint64_t)(beta)*A_2[(i-n_1)*p + k] * B[k*m + j];
                } else if (i < n_1 + n_2 + n_3){
                    ref += (uint64_t)(beta)*A_3[(i-n_1-n_2)*p + k] * B[k*m + j];
                } else {
                    ref += (uint64_t)(beta)*A_4[(i-n_1-n_2-n_3)*p + k] * B[k*m + j];
                }
            }

            if (D_unpacked[i*m + j] != ref) {
                printf("failed\n");
                exit(1);
            }
        }

    }

    printf("passed\n");

    free(A_1);
    free(A_2);
    free(A_3);
    free(A_4);
    free(B);
    free(C);
    free(D_unpacked);
}


// execute LUT mmm per config
template<int n, int m, int p, int num_bits, int n_u, int m_u, int p_u, int N_U, int M_U, int P_U>
void do_lut_mmm() {

    using mat_type = uint8_t;
    using align_type = uint32_t;

    mat_type* A_init = (mat_type*)malloc(n*p*sizeof(mat_type));
    mat_type* B_init = (mat_type*)malloc(p*m*sizeof(mat_type));
    mat_type* C_init = (mat_type*)malloc(n*m*sizeof(mat_type));

    fill<mat_type>(A_init, n, p, num_bits);
    fill<mat_type>(B_init, p, m, num_bits);
    fill<mat_type>(C_init, n, m, num_bits);

    APIB_Mat<num_bits, n, p, mat_type> mat_A(A_init);
    APIB_Mat<num_bits, p, m, mat_type> mat_B(B_init);
    APIB_Mat<num_bits, n, m, mat_type> mat_C(C_init);

    align_type* lut = generate_lut<align_type>(n_u, m_u, p_u, num_bits);

    test_lut(lut, n_u, m_u, p_u, num_bits, 10);

    auto mat_D = lut_mmm<align_type>(mat_A, mat_B, mat_C, N_U, M_U, P_U, n_u, m_u, p_u, lut);
    
    test_mmm<align_type>(mat_D, A_init, B_init, C_init, p);

    free(A_init);
    free(B_init);
    free(C_init);
    free(lut);
}

// mmm tests for generic and slices matrices
// makes sure all cleanup code gets activated
void test_mmm() {
    printf("running mmm tests\n");
    printf("8->32\n");
    do_mmm<15, 701, 200, 2, 2, 2, 2, 2>();
    do_mmm<15, 701, 201, 3, 4, 2, 2, 8>();
    do_mmm<15, 701, 202, 3, 3, 3, 3, 3>();
    do_mmm<15, 701, 203, 3, 4, 2, 3, 6>();
    do_mmm<601, 701, 200, 2, 2, 2, 2, 2>();
    do_mmm<601, 701, 201, 3, 4, 2, 2, 8>();
    do_mmm<601, 701, 202, 3, 3, 3, 3, 3>();
    do_mmm<601, 701, 203, 3, 4, 2, 3, 6>();
    printf("16->32\n");
    do_mmm<1, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<1, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<2, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<2, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<3, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<3, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<4, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<4, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<5, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<5, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<6, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<6, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<7, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<7, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<8, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<8, 701, 201, 6, 7, 2, 17, 11>();
    do_mmm<601, 701, 242, 2, 2, 9, 10, 11>();
    do_mmm<601, 701, 201, 6, 7, 2, 17, 11>();
    printf("64->64 IFMA\n");
    do_mmm<1, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<2, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<3, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<4, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<5, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<6, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<7, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<8, 513, 201, 6, 7, 33, 14, 33>();
    do_mmm<569, 513, 201, 6, 7, 33, 14, 33>();
    printf("64->64\n");
    do_mmm<1, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<2, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<3, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<4, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<5, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<6, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<7, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<8, 513, 201, 6, 7, 33, 18, 33>();
    do_mmm<569, 513, 201, 6, 7, 33, 18, 33>();
    
    // makes sure there are identical slices present
    // also different base types of the slices and different cleanup schedules
    printf("running slice mmm tests\n");
    printf("2 slices\n");
    do_slice_mmm_2<500, 500, 1000, 1000, 2, 2, 2, 3, 4, 4>();
    do_slice_mmm_2<500, 500, 1000, 1000, 2, 2, 10, 3, 4, 4>();
    do_slice_mmm_2<500, 500, 1000, 1000, 2, 2, 33, 33, 4, 4>();
    do_slice_mmm_2<500, 500, 1000, 1000, 2, 2, 33, 33, 19, 33>();
    do_slice_mmm_2<500, 500, 1000, 1000, 2, 2, 10, 33, 19, 33>();
    printf("3 slices\n");
    do_slice_mmm_3<500, 100, 900, 1000, 1000, 2, 2, 2, 3, 4, 4, 4>();
    do_slice_mmm_3<500, 100, 900, 1000, 1000, 2, 2, 10, 3, 5, 4, 4>();
    do_slice_mmm_3<500, 100, 900, 1000, 1000, 2, 2, 33, 33, 30, 4, 4>();
    do_slice_mmm_3<500, 100, 900, 1000, 1000, 2, 2, 33, 33, 32, 19, 33>();
    printf("4 slices\n");
    do_slice_mmm_4<500, 100, 900, 200, 1000, 1000, 2, 2, 2, 3, 4, 5, 4, 4>();
    do_slice_mmm_4<500, 100, 900, 200, 1000, 1000, 2, 2, 10, 3, 4, 5, 4, 4>();
    do_slice_mmm_4<500, 100, 900, 200, 1000, 1000, 2, 2, 33, 33, 32, 30, 4, 4>();
    do_slice_mmm_4<500, 100, 900, 200, 1000, 1000, 2, 2, 33, 33, 30, 32, 19, 33>();

    // check lut mmm
    printf("running lut mmm tests\n");
    do_lut_mmm<384, 384, 384, 1, 4, 4, 2, 128, 128, 128>();
    do_lut_mmm<1024, 1024, 1024, 1, 4, 4, 2, 128, 128, 128>();
}

int main() {
    srand(0);
    test_packing();
    test_unpacking();
    test_scale_unpacking();
    test_mmm();

    printf("passed all tests\n");
}
