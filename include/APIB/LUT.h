/**
 * @file LUT.h
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains functions to build LUT and execute mmm with LUT
 * 
 */

#pragma once

#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <cstdint>
#include <immintrin.h>
#include "APIB.h"
#include "APIB_helper.h"

// utility to compute indices
inline uint64_t compute_address_mmm(uint8_t* a, uint8_t* b, int n_u, int m_u, int p_u, int num_bits, int n, int m, int p) {
    uint64_t res = 0;
    int bit_offset = 0; 
    for (int col = 0; col < m_u; col++) {
        for (int row = 0; row < p_u; row++) {
            res |= (b[row*m + col] << bit_offset);
            bit_offset += num_bits;
        }
    }

    for (int col = 0; col < p_u; col++) {
        for (int row = 0; row < n_u; row++) {
            res |= (a[row*p + col] << bit_offset);
            bit_offset += num_bits;
        }
    }

    return res << ((int)log2(n_u*m_u));
}

// computes the index of a give A-B pair in the table when building and testing
inline uint64_t compute_address(uint8_t* a, uint8_t* b, int n, int m, int p, int num_bits) {
    return compute_address_mmm(a, b, n, m, p, num_bits, n, m, p);
}

// compute all possible B matrices after A
template <typename align_type>
void rec_runner_b(align_type* table, int limit, uint8_t* a, uint8_t* b, int el, int n, int m, int p, int num_bits) {
    // done with filling a and b
    if (el >= p*m) {
        align_type* start = table + compute_address(a, b, n, m, p, num_bits);
        // mmm
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                start[i*m + j] = 0;
                for (int k = 0; k < p; k++) {
                    start[i*m + j] += a[i*p + k] * b[k*m + j];
                }
            }
        }
        // done
        return;
    }

    for (int i = 0; i < limit; i++) {
        b[el] = i;
        rec_runner_b<align_type>(table, limit, a, b, el+1, n, m, p, num_bits);
    }
}

// compute all possible A matrices
template <typename align_type>
void rec_runner_a(align_type* table, int limit, uint8_t* a, uint8_t* b, int el, int n, int m, int p, int num_bits) {
    // done with filling a
    if (el >= n*p) {
        rec_runner_b<align_type>(table, limit, a, b, 0, n, m, p, num_bits);
        return;
    }

    for (int i = 0; i < limit; i++) {
        a[el] = i;
        rec_runner_a<align_type>(table, limit, a, b, el+1, n, m, p, num_bits);
    }
}

// utility to build LUT
template <typename align_type>
void create_lut(align_type* table, int num_bits, int n, int m, int p) {
    int limit = 1 << num_bits;
    uint8_t a[n*p];
    uint8_t b[p*m];

    rec_runner_a<align_type>(table, limit, a, b, 0, n, m, p, num_bits);
}

// test building of LUT
template <typename align_type>
void test_lut(align_type* lut, int n_u, int m_u, int p_u, int num_bits, int num_test_runs) {
    printf("lut build test: ");
    uint8_t a[n_u*p_u];
    uint8_t b[p_u*m_u];
    uint8_t c[p_u*m_u];

    for (int run = 0; run < num_test_runs; run++) {

        for (int i = 0; i < n_u*p_u; i++) {
            a[i] = rand() & ((1 << num_bits) - 1);
        }

        for (int i = 0; i < p_u*m_u; i++) {
            b[i] = rand() & ((1 << num_bits) - 1);
        }

        int start = compute_address(a, b, n_u, m_u, p_u, num_bits);

        for (int i = 0; i < n_u; i++) {
            for (int j = 0; j < m_u; j++) {
                align_type res = 0;
                for (int k = 0; k < p_u; k++) {
                    res += a[i*p_u + k] * b[k*m_u + j];
                }

                if (lut[start + i*m_u + j] != res) {
                    printf("fail!\n");
                    return;
                }
            }
        }

    }

    printf("passed\n");
}

// builds LUT at runtime
template<typename align_type>
align_type* generate_lut(const int n_u, const int m_u, const int p_u, const int num_bits) {
    uint64_t num_combs = (1L << (num_bits*p_u*(n_u+m_u)));
    uint64_t num_els = n_u*m_u*num_combs;

    align_type* lut = (align_type*)aligned_alloc(4096, sizeof(align_type) * num_els);
    create_lut<align_type>(lut, num_bits, n_u, m_u, p_u);
    return lut;
}

// execute mmm using LUT
template<typename align_type, int num_bits_A, int num_bits_B, int num_bits_C, int n, int m, int p, typename type_A, typename type_B, typename type_C>
auto lut_mmm(const APIB_Mat<num_bits_A, n, p, type_A>& mat_A, 
        const APIB_Mat<num_bits_B, p, m, type_B>& mat_B, 
        const APIB_Mat<num_bits_C, n, m, type_C>& mat_C,
        const int N_U, const int M_U, const int P_U,
        const int n_u, const int m_u, const int p_u,
        align_type* lut) {

    const int res_bits = Min_Bits_Res<1, 1, num_bits_A, num_bits_B, num_bits_C, p>::value;

    APIB_Mat<res_bits, n, m, align_type> mat_D; 
    mat_D.reset();

    // hold LUT indices
    uint32_t lut_addr_0[16];
    uint32_t lut_addr_1[16];
    uint32_t lut_addr_2[16];
    uint32_t lut_addr_3[16];

    align_type* c_buf = (align_type*)aligned_alloc(4096, N_U*M_U*sizeof(align_type));

    const uint64_t mask_A = 0b1111;
    auto mask_A_vec = _mm512_set1_epi32(mask_A);
    const int mask_B = 0b11;
    auto mask_B_vec = _mm512_set1_epi32(mask_B);
    auto shift_B_vec_0 = _mm512_set_epi32(  6, 6, 6, 6,
            4, 4, 4, 4,
            2, 2, 2, 2,
            0, 0, 0, 0);

    auto shift_B_vec_1 = _mm512_set_epi32(  6, 4, 2, 0,
            6, 4, 2, 0,
            6, 4, 2, 0,
            6, 4, 2, 0);
    for (int i = 0; i < n; i+=N_U) {
        for (int j = 0; j < m; j+=M_U) {
            unpack<align_type, typename get_base_type<1>::type>(c_buf, mat_C.packed_data, mat_C.get_loc(i, j), N_U, M_U, M_U, num_bits_C);
            for (int k = 0; k < p; k+=P_U) {
                Packed_Loc loc_A = mat_A.get_loc(i, k);
                Packed_Loc loc_B = mat_B.get_loc(k, j);
                uint32_t* mat_A_int = mat_A.packed_data + loc_A.packed_num;
                uint32_t* mat_B_int = mat_B.packed_data + loc_B.packed_num;

                int shift_A = 0;
                for (int i1 = 0; i1 < N_U; i1+=n_u) {
                    for (int j1 = 0; j1 < M_U; j1+=(m_u << 2)) {
                        int addr_count = 16;
                        auto sum_0 = _mm512_setzero_si512();
                        auto sum_1 = _mm512_setzero_si512();
                        auto sum_2 = _mm512_setzero_si512();
                        auto sum_3 = _mm512_setzero_si512();
                        for (int k1 = 0; k1 < P_U; k1+=p_u) {
                            // check if addresses have to be computed
                            if (addr_count >= 16) {
                                addr_count = 0;

                                // A 
                                auto a_vec_0 = _mm512_loadu_si512((void*)(mat_A_int + k1));
                                auto a_vec_1 = _mm512_loadu_si512((void*)(mat_A_int + k1 + 16));

                                auto a_vec_0_shift = _mm512_srli_epi32(a_vec_0, shift_A);
                                auto a_vec_1_shift = _mm512_srli_epi32(a_vec_1, shift_A);

                                auto a_vec_0_mask = _mm512_and_si512(mask_A_vec, a_vec_0_shift);
                                auto a_vec_1_mask = _mm512_and_si512(mask_A_vec, a_vec_1_shift);

                                auto a_vec_0_0 = _mm512_slli_epi32(a_vec_0_mask, 12);
                                auto a_vec_0_1 = _mm512_srli_epi64(a_vec_0_mask, 16);

                                auto a_vec_1_0 = _mm512_slli_epi32(a_vec_1_mask, 12);
                                auto a_vec_1_1 = _mm512_srli_epi64(a_vec_1_mask, 16);

                                auto a_vec_comb_0 = _mm512_or_si512(a_vec_0_0, a_vec_0_1);
                                auto a_vec_comb_1 = _mm512_or_si512(a_vec_1_0, a_vec_1_1);

                                auto a_vec_dense_0 = _mm512_cvtepi64_epi32(a_vec_comb_0);
                                auto a_vec_dense_1 = _mm512_cvtepi64_epi32(a_vec_comb_1);

                                auto a_vec = _mm512_inserti64x4(_mm512_castsi256_si512(a_vec_dense_0), a_vec_dense_1, 1);

                                uint32_t* mat_B_int_block = mat_B_int + (k1 >> 5)*loc_B.eff_cols + j1;

                                // B0
                                auto b_vec_128_0_0 = _mm_loadu_epi32((void*)(mat_B_int_block));
                                auto b_vec_128_1_0 = _mm_srli_epi32(b_vec_128_0_0, 8);
                                auto b_vec_128_2_0 = _mm_srli_epi32(b_vec_128_0_0, 16);
                                auto b_vec_128_3_0 = _mm_srli_epi32(b_vec_128_0_0, 24);

                                auto b_vec_broad_0_0 = _mm512_broadcast_i32x4(b_vec_128_0_0);
                                auto b_vec_broad_1_0 = _mm512_broadcast_i32x4(b_vec_128_1_0);
                                auto b_vec_broad_2_0 = _mm512_broadcast_i32x4(b_vec_128_2_0);
                                auto b_vec_broad_3_0 = _mm512_broadcast_i32x4(b_vec_128_3_0);

                                auto b_vec_0_shifted_0 = _mm512_srlv_epi32(b_vec_broad_0_0, shift_B_vec_0);
                                auto b_vec_1_shifted_0 = _mm512_srlv_epi32(b_vec_broad_1_0, shift_B_vec_0);
                                auto b_vec_2_shifted_0 = _mm512_srlv_epi32(b_vec_broad_2_0, shift_B_vec_0);
                                auto b_vec_3_shifted_0 = _mm512_srlv_epi32(b_vec_broad_3_0, shift_B_vec_0);

                                auto b_vec_0_0 = _mm512_and_si512(mask_B_vec, b_vec_0_shifted_0);
                                auto b_vec_1_0 = _mm512_and_si512(mask_B_vec, b_vec_1_shifted_0);
                                auto b_vec_2_0 = _mm512_and_si512(mask_B_vec, b_vec_2_shifted_0);
                                auto b_vec_3_0 = _mm512_and_si512(mask_B_vec, b_vec_3_shifted_0);

                                auto b_vec_0_in_pos_0 = _mm512_sllv_epi32(b_vec_0_0, shift_B_vec_1);
                                auto b_vec_1_in_pos_0 = _mm512_sllv_epi32(b_vec_1_0, shift_B_vec_1);
                                auto b_vec_2_in_pos_0 = _mm512_sllv_epi32(b_vec_2_0, shift_B_vec_1);
                                auto b_vec_3_in_pos_0 = _mm512_sllv_epi32(b_vec_3_0, shift_B_vec_1);

                                auto b_vec_dense_0_0_0 = b_vec_0_in_pos_0;
                                auto b_vec_dense_0_1_0 = _mm512_bsrli_epi128(b_vec_0_in_pos_0, 4);
                                auto b_vec_dense_0_2_0 = _mm512_bsrli_epi128(b_vec_0_in_pos_0, 8);
                                auto b_vec_dense_0_3_0 = _mm512_bsrli_epi128(b_vec_0_in_pos_0, 12);

                                auto b_vec_dense_1_0_0 = b_vec_1_in_pos_0;
                                auto b_vec_dense_1_1_0 = _mm512_bsrli_epi128(b_vec_1_in_pos_0, 4);
                                auto b_vec_dense_1_2_0 = _mm512_bsrli_epi128(b_vec_1_in_pos_0, 8);
                                auto b_vec_dense_1_3_0 = _mm512_bsrli_epi128(b_vec_1_in_pos_0, 12);

                                auto b_vec_dense_2_0_0 = b_vec_2_in_pos_0;
                                auto b_vec_dense_2_1_0 = _mm512_bsrli_epi128(b_vec_2_in_pos_0, 4);
                                auto b_vec_dense_2_2_0 = _mm512_bsrli_epi128(b_vec_2_in_pos_0, 8);
                                auto b_vec_dense_2_3_0 = _mm512_bsrli_epi128(b_vec_2_in_pos_0, 12);

                                auto b_vec_dense_3_0_0 = b_vec_3_in_pos_0;
                                auto b_vec_dense_3_1_0 = _mm512_bsrli_epi128(b_vec_3_in_pos_0, 4);
                                auto b_vec_dense_3_2_0 = _mm512_bsrli_epi128(b_vec_3_in_pos_0, 8);
                                auto b_vec_dense_3_3_0 = _mm512_bsrli_epi128(b_vec_3_in_pos_0, 12);

                                auto b_vec_comb_0_0_1_0 = _mm512_or_si512(b_vec_dense_0_0_0, b_vec_dense_0_1_0);
                                auto b_vec_comb_0_2_3_0 = _mm512_or_si512(b_vec_dense_0_2_0, b_vec_dense_0_3_0);

                                auto b_vec_comb_1_0_1_0 = _mm512_or_si512(b_vec_dense_1_0_0, b_vec_dense_1_1_0);
                                auto b_vec_comb_1_2_3_0 = _mm512_or_si512(b_vec_dense_1_2_0, b_vec_dense_1_3_0);

                                auto b_vec_comb_2_0_1_0 = _mm512_or_si512(b_vec_dense_2_0_0, b_vec_dense_2_1_0);
                                auto b_vec_comb_2_2_3_0 = _mm512_or_si512(b_vec_dense_2_2_0, b_vec_dense_2_3_0);

                                auto b_vec_comb_3_0_1_0 = _mm512_or_si512(b_vec_dense_3_0_0, b_vec_dense_3_1_0);
                                auto b_vec_comb_3_2_3_0 = _mm512_or_si512(b_vec_dense_3_2_0, b_vec_dense_3_3_0);

                                auto b_vec_comb_0_0 = _mm512_or_si512(b_vec_comb_0_0_1_0, b_vec_comb_0_2_3_0);
                                auto b_vec_comb_1_0 = _mm512_or_si512(b_vec_comb_1_0_1_0, b_vec_comb_1_2_3_0);
                                auto b_vec_comb_2_0 = _mm512_or_si512(b_vec_comb_2_0_1_0, b_vec_comb_2_2_3_0);
                                auto b_vec_comb_3_0 = _mm512_or_si512(b_vec_comb_3_0_1_0, b_vec_comb_3_2_3_0);

                                auto b_vec_comb_0_dense_0 = _mm512_cvtepi64_epi32(b_vec_comb_0_0);
                                auto b_vec_comb_1_dense_0 = _mm512_cvtepi64_epi32(b_vec_comb_1_0);
                                auto b_vec_comb_2_dense_0 = _mm512_cvtepi64_epi32(b_vec_comb_2_0);
                                auto b_vec_comb_3_dense_0 = _mm512_cvtepi64_epi32(b_vec_comb_3_0);

                                auto b_vec_comb_0_1_dense_0 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_dense_0), b_vec_comb_1_dense_0, 1);
                                auto b_vec_comb_2_3_dense_0 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_2_dense_0), b_vec_comb_3_dense_0, 1);

                                auto b_vec_comb_0_1_dense_dense_0 = _mm512_cvtepi64_epi32(b_vec_comb_0_1_dense_0);
                                auto b_vec_comb_2_3_dense_dense_0 = _mm512_cvtepi64_epi32(b_vec_comb_2_3_dense_0);

                                auto b_vec_dense_0 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_1_dense_dense_0), b_vec_comb_2_3_dense_dense_0, 1);

                                auto b_vec_0 = _mm512_slli_epi32(b_vec_dense_0, 4);

                                auto addr_vec_0 = _mm512_or_si512(a_vec, b_vec_0);

                                _mm512_storeu_si512((void*)lut_addr_0, addr_vec_0);

                                // B1
                                auto b_vec_128_0_1 = _mm_loadu_epi32((void*)(mat_B_int_block + m_u));
                                auto b_vec_128_1_1 = _mm_srli_epi32(b_vec_128_0_1, 8);
                                auto b_vec_128_2_1 = _mm_srli_epi32(b_vec_128_0_1, 16);
                                auto b_vec_128_3_1 = _mm_srli_epi32(b_vec_128_0_1, 24);

                                auto b_vec_broad_0_1 = _mm512_broadcast_i32x4(b_vec_128_0_1);
                                auto b_vec_broad_1_1 = _mm512_broadcast_i32x4(b_vec_128_1_1);
                                auto b_vec_broad_2_1 = _mm512_broadcast_i32x4(b_vec_128_2_1);
                                auto b_vec_broad_3_1 = _mm512_broadcast_i32x4(b_vec_128_3_1);

                                auto b_vec_0_shifted_1 = _mm512_srlv_epi32(b_vec_broad_0_1, shift_B_vec_0);
                                auto b_vec_1_shifted_1 = _mm512_srlv_epi32(b_vec_broad_1_1, shift_B_vec_0);
                                auto b_vec_2_shifted_1 = _mm512_srlv_epi32(b_vec_broad_2_1, shift_B_vec_0);
                                auto b_vec_3_shifted_1 = _mm512_srlv_epi32(b_vec_broad_3_1, shift_B_vec_0);

                                auto b_vec_0_1 = _mm512_and_si512(mask_B_vec, b_vec_0_shifted_1);
                                auto b_vec_1_1 = _mm512_and_si512(mask_B_vec, b_vec_1_shifted_1);
                                auto b_vec_2_1 = _mm512_and_si512(mask_B_vec, b_vec_2_shifted_1);
                                auto b_vec_3_1 = _mm512_and_si512(mask_B_vec, b_vec_3_shifted_1);

                                auto b_vec_0_in_pos_1 = _mm512_sllv_epi32(b_vec_0_1, shift_B_vec_1);
                                auto b_vec_1_in_pos_1 = _mm512_sllv_epi32(b_vec_1_1, shift_B_vec_1);
                                auto b_vec_2_in_pos_1 = _mm512_sllv_epi32(b_vec_2_1, shift_B_vec_1);
                                auto b_vec_3_in_pos_1 = _mm512_sllv_epi32(b_vec_3_1, shift_B_vec_1);

                                auto b_vec_dense_0_0_1 = b_vec_0_in_pos_1;
                                auto b_vec_dense_0_1_1 = _mm512_bsrli_epi128(b_vec_0_in_pos_1, 4);
                                auto b_vec_dense_0_2_1 = _mm512_bsrli_epi128(b_vec_0_in_pos_1, 8);
                                auto b_vec_dense_0_3_1 = _mm512_bsrli_epi128(b_vec_0_in_pos_1, 12);

                                auto b_vec_dense_1_0_1 = b_vec_1_in_pos_1;
                                auto b_vec_dense_1_1_1 = _mm512_bsrli_epi128(b_vec_1_in_pos_1, 4);
                                auto b_vec_dense_1_2_1 = _mm512_bsrli_epi128(b_vec_1_in_pos_1, 8);
                                auto b_vec_dense_1_3_1 = _mm512_bsrli_epi128(b_vec_1_in_pos_1, 12);

                                auto b_vec_dense_2_0_1 = b_vec_2_in_pos_1;
                                auto b_vec_dense_2_1_1 = _mm512_bsrli_epi128(b_vec_2_in_pos_1, 4);
                                auto b_vec_dense_2_2_1 = _mm512_bsrli_epi128(b_vec_2_in_pos_1, 8);
                                auto b_vec_dense_2_3_1 = _mm512_bsrli_epi128(b_vec_2_in_pos_1, 12);

                                auto b_vec_dense_3_0_1 = b_vec_3_in_pos_1;
                                auto b_vec_dense_3_1_1 = _mm512_bsrli_epi128(b_vec_3_in_pos_1, 4);
                                auto b_vec_dense_3_2_1 = _mm512_bsrli_epi128(b_vec_3_in_pos_1, 8);
                                auto b_vec_dense_3_3_1 = _mm512_bsrli_epi128(b_vec_3_in_pos_1, 12);

                                auto b_vec_comb_0_0_1_1 = _mm512_or_si512(b_vec_dense_0_0_1, b_vec_dense_0_1_1);
                                auto b_vec_comb_0_2_3_1 = _mm512_or_si512(b_vec_dense_0_2_1, b_vec_dense_0_3_1);

                                auto b_vec_comb_1_0_1_1 = _mm512_or_si512(b_vec_dense_1_0_1, b_vec_dense_1_1_1);
                                auto b_vec_comb_1_2_3_1 = _mm512_or_si512(b_vec_dense_1_2_1, b_vec_dense_1_3_1);

                                auto b_vec_comb_2_0_1_1 = _mm512_or_si512(b_vec_dense_2_0_1, b_vec_dense_2_1_1);
                                auto b_vec_comb_2_2_3_1 = _mm512_or_si512(b_vec_dense_2_2_1, b_vec_dense_2_3_1);

                                auto b_vec_comb_3_0_1_1 = _mm512_or_si512(b_vec_dense_3_0_1, b_vec_dense_3_1_1);
                                auto b_vec_comb_3_2_3_1 = _mm512_or_si512(b_vec_dense_3_2_1, b_vec_dense_3_3_1);

                                auto b_vec_comb_0_1 = _mm512_or_si512(b_vec_comb_0_0_1_1, b_vec_comb_0_2_3_1);
                                auto b_vec_comb_1_1 = _mm512_or_si512(b_vec_comb_1_0_1_1, b_vec_comb_1_2_3_1);
                                auto b_vec_comb_2_1 = _mm512_or_si512(b_vec_comb_2_0_1_1, b_vec_comb_2_2_3_1);
                                auto b_vec_comb_3_1 = _mm512_or_si512(b_vec_comb_3_0_1_1, b_vec_comb_3_2_3_1);

                                auto b_vec_comb_0_dense_1 = _mm512_cvtepi64_epi32(b_vec_comb_0_1);
                                auto b_vec_comb_1_dense_1 = _mm512_cvtepi64_epi32(b_vec_comb_1_1);
                                auto b_vec_comb_2_dense_1 = _mm512_cvtepi64_epi32(b_vec_comb_2_1);
                                auto b_vec_comb_3_dense_1 = _mm512_cvtepi64_epi32(b_vec_comb_3_1);

                                auto b_vec_comb_0_1_dense_1 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_dense_1), b_vec_comb_1_dense_1, 1);
                                auto b_vec_comb_2_3_dense_1 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_2_dense_1), b_vec_comb_3_dense_1, 1);

                                auto b_vec_comb_0_1_dense_dense_1 = _mm512_cvtepi64_epi32(b_vec_comb_0_1_dense_1);
                                auto b_vec_comb_2_3_dense_dense_1 = _mm512_cvtepi64_epi32(b_vec_comb_2_3_dense_1);

                                auto b_vec_dense_1 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_1_dense_dense_1), b_vec_comb_2_3_dense_dense_1, 1);

                                auto b_vec_1 = _mm512_slli_epi32(b_vec_dense_1, 4);

                                auto addr_vec_1 = _mm512_or_si512(a_vec, b_vec_1);

                                _mm512_storeu_si512((void*)lut_addr_1, addr_vec_1);

                                // B2
                                auto b_vec_128_0_2 = _mm_loadu_epi32((void*)(mat_B_int_block + 2*m_u));
                                auto b_vec_128_1_2 = _mm_srli_epi32(b_vec_128_0_2, 8);
                                auto b_vec_128_2_2 = _mm_srli_epi32(b_vec_128_0_2, 16);
                                auto b_vec_128_3_2 = _mm_srli_epi32(b_vec_128_0_2, 24);

                                auto b_vec_broad_0_2 = _mm512_broadcast_i32x4(b_vec_128_0_2);
                                auto b_vec_broad_1_2 = _mm512_broadcast_i32x4(b_vec_128_1_2);
                                auto b_vec_broad_2_2 = _mm512_broadcast_i32x4(b_vec_128_2_2);
                                auto b_vec_broad_3_2 = _mm512_broadcast_i32x4(b_vec_128_3_2);

                                auto b_vec_0_shifted_2 = _mm512_srlv_epi32(b_vec_broad_0_2, shift_B_vec_0);
                                auto b_vec_1_shifted_2 = _mm512_srlv_epi32(b_vec_broad_1_2, shift_B_vec_0);
                                auto b_vec_2_shifted_2 = _mm512_srlv_epi32(b_vec_broad_2_2, shift_B_vec_0);
                                auto b_vec_3_shifted_2 = _mm512_srlv_epi32(b_vec_broad_3_2, shift_B_vec_0);

                                auto b_vec_0_2 = _mm512_and_si512(mask_B_vec, b_vec_0_shifted_2);
                                auto b_vec_1_2 = _mm512_and_si512(mask_B_vec, b_vec_1_shifted_2);
                                auto b_vec_2_2 = _mm512_and_si512(mask_B_vec, b_vec_2_shifted_2);
                                auto b_vec_3_2 = _mm512_and_si512(mask_B_vec, b_vec_3_shifted_2);

                                auto b_vec_0_in_pos_2 = _mm512_sllv_epi32(b_vec_0_2, shift_B_vec_1);
                                auto b_vec_1_in_pos_2 = _mm512_sllv_epi32(b_vec_1_2, shift_B_vec_1);
                                auto b_vec_2_in_pos_2 = _mm512_sllv_epi32(b_vec_2_2, shift_B_vec_1);
                                auto b_vec_3_in_pos_2 = _mm512_sllv_epi32(b_vec_3_2, shift_B_vec_1);

                                auto b_vec_dense_0_0_2 = b_vec_0_in_pos_2;
                                auto b_vec_dense_0_1_2 = _mm512_bsrli_epi128(b_vec_0_in_pos_2, 4);
                                auto b_vec_dense_0_2_2 = _mm512_bsrli_epi128(b_vec_0_in_pos_2, 8);
                                auto b_vec_dense_0_3_2 = _mm512_bsrli_epi128(b_vec_0_in_pos_2, 12);

                                auto b_vec_dense_1_0_2 = b_vec_1_in_pos_2;
                                auto b_vec_dense_1_1_2 = _mm512_bsrli_epi128(b_vec_1_in_pos_2, 4);
                                auto b_vec_dense_1_2_2 = _mm512_bsrli_epi128(b_vec_1_in_pos_2, 8);
                                auto b_vec_dense_1_3_2 = _mm512_bsrli_epi128(b_vec_1_in_pos_2, 12);

                                auto b_vec_dense_2_0_2 = b_vec_2_in_pos_2;
                                auto b_vec_dense_2_1_2 = _mm512_bsrli_epi128(b_vec_2_in_pos_2, 4);
                                auto b_vec_dense_2_2_2 = _mm512_bsrli_epi128(b_vec_2_in_pos_2, 8);
                                auto b_vec_dense_2_3_2 = _mm512_bsrli_epi128(b_vec_2_in_pos_2, 12);

                                auto b_vec_dense_3_0_2 = b_vec_3_in_pos_2;
                                auto b_vec_dense_3_1_2 = _mm512_bsrli_epi128(b_vec_3_in_pos_2, 4);
                                auto b_vec_dense_3_2_2 = _mm512_bsrli_epi128(b_vec_3_in_pos_2, 8);
                                auto b_vec_dense_3_3_2 = _mm512_bsrli_epi128(b_vec_3_in_pos_2, 12);

                                auto b_vec_comb_0_0_1_2 = _mm512_or_si512(b_vec_dense_0_0_2, b_vec_dense_0_1_2);
                                auto b_vec_comb_0_2_3_2 = _mm512_or_si512(b_vec_dense_0_2_2, b_vec_dense_0_3_2);

                                auto b_vec_comb_1_0_1_2 = _mm512_or_si512(b_vec_dense_1_0_2, b_vec_dense_1_1_2);
                                auto b_vec_comb_1_2_3_2 = _mm512_or_si512(b_vec_dense_1_2_2, b_vec_dense_1_3_2);

                                auto b_vec_comb_2_0_1_2 = _mm512_or_si512(b_vec_dense_2_0_2, b_vec_dense_2_1_2);
                                auto b_vec_comb_2_2_3_2 = _mm512_or_si512(b_vec_dense_2_2_2, b_vec_dense_2_3_2);

                                auto b_vec_comb_3_0_1_2 = _mm512_or_si512(b_vec_dense_3_0_2, b_vec_dense_3_1_2);
                                auto b_vec_comb_3_2_3_2 = _mm512_or_si512(b_vec_dense_3_2_2, b_vec_dense_3_3_2);

                                auto b_vec_comb_0_2 = _mm512_or_si512(b_vec_comb_0_0_1_2, b_vec_comb_0_2_3_2);
                                auto b_vec_comb_1_2 = _mm512_or_si512(b_vec_comb_1_0_1_2, b_vec_comb_1_2_3_2);
                                auto b_vec_comb_2_2 = _mm512_or_si512(b_vec_comb_2_0_1_2, b_vec_comb_2_2_3_2);
                                auto b_vec_comb_3_2 = _mm512_or_si512(b_vec_comb_3_0_1_2, b_vec_comb_3_2_3_2);

                                auto b_vec_comb_0_dense_2 = _mm512_cvtepi64_epi32(b_vec_comb_0_2);
                                auto b_vec_comb_1_dense_2 = _mm512_cvtepi64_epi32(b_vec_comb_1_2);
                                auto b_vec_comb_2_dense_2 = _mm512_cvtepi64_epi32(b_vec_comb_2_2);
                                auto b_vec_comb_3_dense_2 = _mm512_cvtepi64_epi32(b_vec_comb_3_2);

                                auto b_vec_comb_0_1_dense_2 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_dense_2), b_vec_comb_1_dense_2, 1);
                                auto b_vec_comb_2_3_dense_2 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_2_dense_2), b_vec_comb_3_dense_2, 1);

                                auto b_vec_comb_0_1_dense_dense_2 = _mm512_cvtepi64_epi32(b_vec_comb_0_1_dense_2);
                                auto b_vec_comb_2_3_dense_dense_2 = _mm512_cvtepi64_epi32(b_vec_comb_2_3_dense_2);

                                auto b_vec_dense_2 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_1_dense_dense_2), b_vec_comb_2_3_dense_dense_2, 1);

                                auto b_vec_2 = _mm512_slli_epi32(b_vec_dense_2, 4);

                                auto addr_vec_2 = _mm512_or_si512(a_vec, b_vec_2);

                                _mm512_storeu_si512((void*)lut_addr_2, addr_vec_2);

                                // B3
                                auto b_vec_128_0_3 = _mm_loadu_epi32((void*)(mat_B_int_block + 3*m_u));
                                auto b_vec_128_1_3 = _mm_srli_epi32(b_vec_128_0_3, 8);
                                auto b_vec_128_2_3 = _mm_srli_epi32(b_vec_128_0_3, 16);
                                auto b_vec_128_3_3 = _mm_srli_epi32(b_vec_128_0_3, 24);

                                auto b_vec_broad_0_3 = _mm512_broadcast_i32x4(b_vec_128_0_3);
                                auto b_vec_broad_1_3 = _mm512_broadcast_i32x4(b_vec_128_1_3);
                                auto b_vec_broad_2_3 = _mm512_broadcast_i32x4(b_vec_128_2_3);
                                auto b_vec_broad_3_3 = _mm512_broadcast_i32x4(b_vec_128_3_3);

                                auto b_vec_0_shifted_3 = _mm512_srlv_epi32(b_vec_broad_0_3, shift_B_vec_0);
                                auto b_vec_1_shifted_3 = _mm512_srlv_epi32(b_vec_broad_1_3, shift_B_vec_0);
                                auto b_vec_2_shifted_3 = _mm512_srlv_epi32(b_vec_broad_2_3, shift_B_vec_0);
                                auto b_vec_3_shifted_3 = _mm512_srlv_epi32(b_vec_broad_3_3, shift_B_vec_0);

                                auto b_vec_0_3 = _mm512_and_si512(mask_B_vec, b_vec_0_shifted_3);
                                auto b_vec_1_3 = _mm512_and_si512(mask_B_vec, b_vec_1_shifted_3);
                                auto b_vec_2_3 = _mm512_and_si512(mask_B_vec, b_vec_2_shifted_3);
                                auto b_vec_3_3 = _mm512_and_si512(mask_B_vec, b_vec_3_shifted_3);

                                auto b_vec_0_in_pos_3 = _mm512_sllv_epi32(b_vec_0_3, shift_B_vec_1);
                                auto b_vec_1_in_pos_3 = _mm512_sllv_epi32(b_vec_1_3, shift_B_vec_1);
                                auto b_vec_2_in_pos_3 = _mm512_sllv_epi32(b_vec_2_3, shift_B_vec_1);
                                auto b_vec_3_in_pos_3 = _mm512_sllv_epi32(b_vec_3_3, shift_B_vec_1);

                                auto b_vec_dense_0_0_3 = b_vec_0_in_pos_3;
                                auto b_vec_dense_0_1_3 = _mm512_bsrli_epi128(b_vec_0_in_pos_3, 4);
                                auto b_vec_dense_0_2_3 = _mm512_bsrli_epi128(b_vec_0_in_pos_3, 8);
                                auto b_vec_dense_0_3_3 = _mm512_bsrli_epi128(b_vec_0_in_pos_3, 12);

                                auto b_vec_dense_1_0_3 = b_vec_1_in_pos_3;
                                auto b_vec_dense_1_1_3 = _mm512_bsrli_epi128(b_vec_1_in_pos_3, 4);
                                auto b_vec_dense_1_2_3 = _mm512_bsrli_epi128(b_vec_1_in_pos_3, 8);
                                auto b_vec_dense_1_3_3 = _mm512_bsrli_epi128(b_vec_1_in_pos_3, 12);

                                auto b_vec_dense_2_0_3 = b_vec_2_in_pos_3;
                                auto b_vec_dense_2_1_3 = _mm512_bsrli_epi128(b_vec_2_in_pos_3, 4);
                                auto b_vec_dense_2_2_3 = _mm512_bsrli_epi128(b_vec_2_in_pos_3, 8);
                                auto b_vec_dense_2_3_3 = _mm512_bsrli_epi128(b_vec_2_in_pos_3, 12);

                                auto b_vec_dense_3_0_3 = b_vec_3_in_pos_3;
                                auto b_vec_dense_3_1_3 = _mm512_bsrli_epi128(b_vec_3_in_pos_3, 4);
                                auto b_vec_dense_3_2_3 = _mm512_bsrli_epi128(b_vec_3_in_pos_3, 8);
                                auto b_vec_dense_3_3_3 = _mm512_bsrli_epi128(b_vec_3_in_pos_3, 12);

                                auto b_vec_comb_0_0_1_3 = _mm512_or_si512(b_vec_dense_0_0_3, b_vec_dense_0_1_3);
                                auto b_vec_comb_0_2_3_3 = _mm512_or_si512(b_vec_dense_0_2_3, b_vec_dense_0_3_3);

                                auto b_vec_comb_1_0_1_3 = _mm512_or_si512(b_vec_dense_1_0_3, b_vec_dense_1_1_3);
                                auto b_vec_comb_1_2_3_3 = _mm512_or_si512(b_vec_dense_1_2_3, b_vec_dense_1_3_3);

                                auto b_vec_comb_2_0_1_3 = _mm512_or_si512(b_vec_dense_2_0_3, b_vec_dense_2_1_3);
                                auto b_vec_comb_2_2_3_3 = _mm512_or_si512(b_vec_dense_2_2_3, b_vec_dense_2_3_3);

                                auto b_vec_comb_3_0_1_3 = _mm512_or_si512(b_vec_dense_3_0_3, b_vec_dense_3_1_3);
                                auto b_vec_comb_3_2_3_3 = _mm512_or_si512(b_vec_dense_3_2_3, b_vec_dense_3_3_3);

                                auto b_vec_comb_0_3 = _mm512_or_si512(b_vec_comb_0_0_1_3, b_vec_comb_0_2_3_3);
                                auto b_vec_comb_1_3 = _mm512_or_si512(b_vec_comb_1_0_1_3, b_vec_comb_1_2_3_3);
                                auto b_vec_comb_2_3 = _mm512_or_si512(b_vec_comb_2_0_1_3, b_vec_comb_2_2_3_3);
                                auto b_vec_comb_3_3 = _mm512_or_si512(b_vec_comb_3_0_1_3, b_vec_comb_3_2_3_3);

                                auto b_vec_comb_0_dense_3 = _mm512_cvtepi64_epi32(b_vec_comb_0_3);
                                auto b_vec_comb_1_dense_3 = _mm512_cvtepi64_epi32(b_vec_comb_1_3);
                                auto b_vec_comb_2_dense_3 = _mm512_cvtepi64_epi32(b_vec_comb_2_3);
                                auto b_vec_comb_3_dense_3 = _mm512_cvtepi64_epi32(b_vec_comb_3_3);

                                auto b_vec_comb_0_1_dense_3 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_dense_3), b_vec_comb_1_dense_3, 1);
                                auto b_vec_comb_2_3_dense_3 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_2_dense_3), b_vec_comb_3_dense_3, 1);

                                auto b_vec_comb_0_1_dense_dense_3 = _mm512_cvtepi64_epi32(b_vec_comb_0_1_dense_3);
                                auto b_vec_comb_2_3_dense_dense_3 = _mm512_cvtepi64_epi32(b_vec_comb_2_3_dense_3);

                                auto b_vec_dense_3 = _mm512_inserti64x4(_mm512_castsi256_si512(b_vec_comb_0_1_dense_dense_3), b_vec_comb_2_3_dense_dense_3, 1);

                                auto b_vec_3 = _mm512_slli_epi32(b_vec_dense_3, 4);

                                auto addr_vec_3 = _mm512_or_si512(a_vec, b_vec_3);

                                _mm512_storeu_si512((void*)lut_addr_3, addr_vec_3);
                            }

                            align_type* lut_mat_0 = lut + lut_addr_0[addr_count];
                            align_type* lut_mat_1 = lut + lut_addr_1[addr_count];
                            align_type* lut_mat_2 = lut + lut_addr_2[addr_count];
                            align_type* lut_mat_3 = lut + lut_addr_3[addr_count];

                            auto mat_vec_0 = _mm512_loadu_si512((void*)lut_mat_0);
                            auto mat_vec_1 = _mm512_loadu_si512((void*)lut_mat_1);
                            auto mat_vec_2 = _mm512_loadu_si512((void*)lut_mat_2);
                            auto mat_vec_3 = _mm512_loadu_si512((void*)lut_mat_3);

                            sum_0 = _mm512_add_epi32(mat_vec_0, sum_0);
                            sum_1 = _mm512_add_epi32(mat_vec_1, sum_1);
                            sum_2 = _mm512_add_epi32(mat_vec_2, sum_2);
                            sum_3 = _mm512_add_epi32(mat_vec_3, sum_3);

                            addr_count++;
                        }

                        // transform
                        const int all_0 = 0b00000000;
                        const int all_1 = 0b01010101;
                        const int all_2 = 0b10101010;
                        const int all_3 = 0b11111111;
                        auto mask_0 = _mm512_int2mask(0b1111000011110000);
                        auto mask_1 = _mm512_int2mask(0b1111111100000000);

                        auto perm_vec_0_0 = sum_0;
                        auto perm_vec_0_1 = _mm512_shuffle_i64x2(sum_0, sum_0, all_1);
                        auto perm_vec_0_2 = _mm512_shuffle_i64x2(sum_0, sum_0, all_2);
                        auto perm_vec_0_3 = _mm512_shuffle_i64x2(sum_0, sum_0, all_3);

                        auto perm_vec_1_0 = _mm512_shuffle_i64x2(sum_1, sum_1, all_0);
                        auto perm_vec_1_1 = sum_1;
                        auto perm_vec_1_2 = _mm512_shuffle_i64x2(sum_1, sum_1, all_2);
                        auto perm_vec_1_3 = _mm512_shuffle_i64x2(sum_1, sum_1, all_3);

                        auto perm_vec_2_0 = _mm512_shuffle_i64x2(sum_2, sum_2, all_0);
                        auto perm_vec_2_1 = _mm512_shuffle_i64x2(sum_2, sum_2, all_1);
                        auto perm_vec_2_2 = sum_2;
                        auto perm_vec_2_3 = _mm512_shuffle_i64x2(sum_2, sum_2, all_3);

                        auto perm_vec_3_0 = _mm512_shuffle_i64x2(sum_3, sum_3, all_0);
                        auto perm_vec_3_1 = _mm512_shuffle_i64x2(sum_3, sum_3, all_1);
                        auto perm_vec_3_2 = _mm512_shuffle_i64x2(sum_3, sum_3, all_2);
                        auto perm_vec_3_3 = sum_3;

                        auto blend_left_0 = _mm512_mask_blend_epi32(mask_0, perm_vec_0_0, perm_vec_1_0);
                        auto blend_left_1 = _mm512_mask_blend_epi32(mask_0, perm_vec_0_1, perm_vec_1_1);
                        auto blend_left_2 = _mm512_mask_blend_epi32(mask_0, perm_vec_0_2, perm_vec_1_2);
                        auto blend_left_3 = _mm512_mask_blend_epi32(mask_0, perm_vec_0_3, perm_vec_1_3);

                        auto blend_right_0 = _mm512_mask_blend_epi32(mask_0, perm_vec_2_0, perm_vec_3_0);
                        auto blend_right_1 = _mm512_mask_blend_epi32(mask_0, perm_vec_2_1, perm_vec_3_1);
                        auto blend_right_2 = _mm512_mask_blend_epi32(mask_0, perm_vec_2_2, perm_vec_3_2);
                        auto blend_right_3 = _mm512_mask_blend_epi32(mask_0, perm_vec_2_3, perm_vec_3_3);

                        auto blend_0 = _mm512_mask_blend_epi32(mask_1, blend_left_0, blend_right_0);
                        auto blend_1 = _mm512_mask_blend_epi32(mask_1, blend_left_1, blend_right_1);
                        auto blend_2 = _mm512_mask_blend_epi32(mask_1, blend_left_2, blend_right_2);
                        auto blend_3 = _mm512_mask_blend_epi32(mask_1, blend_left_3, blend_right_3);

                        auto c_vec_0 = _mm512_loadu_si512((void*)(c_buf + i1*M_U + j1));
                        auto c_vec_1 = _mm512_loadu_si512((void*)(c_buf + (i1+1)*M_U + j1));
                        auto c_vec_2 = _mm512_loadu_si512((void*)(c_buf + (i1+2)*M_U + j1));
                        auto c_vec_3 = _mm512_loadu_si512((void*)(c_buf + (i1+3)*M_U + j1));

                        auto res_vec_0 = _mm512_add_epi32(blend_0, c_vec_0);
                        auto res_vec_1 = _mm512_add_epi32(blend_1, c_vec_1);
                        auto res_vec_2 = _mm512_add_epi32(blend_2, c_vec_2);
                        auto res_vec_3 = _mm512_add_epi32(blend_3, c_vec_3);

                        _mm512_storeu_si512((void*)(c_buf + i1*M_U + j1), res_vec_0);
                        _mm512_storeu_si512((void*)(c_buf + (i1+1)*M_U + j1), res_vec_1);
                        _mm512_storeu_si512((void*)(c_buf + (i1+2)*M_U + j1), res_vec_2);
                        _mm512_storeu_si512((void*)(c_buf + (i1+3)*M_U + j1), res_vec_3);
                    }

                    // adjust shift
                    shift_A += 4;
                    if (shift_A == 32) {
                        mat_A_int += loc_A.eff_cols;
                        shift_A = 0;
                    }
                }
            }
            pack<align_type, typename get_base_type<res_bits>::type>(c_buf, mat_D.packed_data, mat_D.get_loc(i, j), N_U, M_U, M_U, res_bits);
        }
    }

    free(c_buf);

    return mat_D;
}

// tests LUT mmm for correctness
template<typename align_type, typename mat_type, int num_bits, int n, int m, typename type>
void test_mmm(APIB_Mat<num_bits, n, m, type>& mat_D, mat_type* A, mat_type* B, mat_type* C, int p) {
    align_type* D_unpacked = (align_type*)malloc(n*m*sizeof(align_type));
    mat_D.template unpack_data<align_type>(D_unpacked, m); 

    printf("lut mmm test: ");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            align_type res = (align_type)C[i*m + j];
            for (int k = 0; k < p; k++) {
                res += (align_type)A[i*p + k] * (align_type)B[k*m + j];
            }

            if (res != D_unpacked[i*m + j]) {
                printf("failed!\n");
                exit(1);
            }
        }
    }
    printf("passed\n");
}
