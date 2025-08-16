/**
 * @file APIB_packing.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains implementations of specializations the pack function
 * 
 */

#include "APIB/APIB_packing.h"
#include "APIB/APIB_helper.h"
#include <cstdint>

/**
 * @brief Packs data from row-major buffer
 * 
 * @tparam buf_type Type of source buffer
 * @tparam base_type Basetype of APIB Mat object containing packed data
 * @param buf Pointer to source buffer
 * @param packed_data Pointer to packed data
 * @param loc Location of offset
 * @param row_num Number of rows to pack
 * @param col_num Number of columns to pack, assumes stays in same block
 * @param buf_cols Number of columns of source buffer
 * @param num_bits Size of elements in bits
 */
template<typename buf_type, typename base_type>
void pack(buf_type* buf, base_type* packed_data, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("generic packing\n");
    int bits_type = 8 * sizeof(base_type);
    for (int row = 0; row < row_num; row++) {
        // go through columns
        for (int col = 0; col < col_num; col++) {
            // store bits
            packed_data[loc.packed_num + col] |= (base_type)(buf[row * buf_cols + col] << loc.bit_offset);
        }

        // adjust current bit position
        loc.bit_offset += num_bits;

        // wrap around
        if (loc.bit_offset >= bits_type) {
            loc.bit_offset &= (bits_type - 1);
            loc.packed_num += loc.eff_cols;

            // initialize for next run
            if (loc.bit_offset != 0) {
                for (int col = 0; col < col_num; col++) {
                    // store remaining bits from last element
                    packed_data[loc.packed_num + col] = (base_type)(buf[row * buf_cols + col] >> (num_bits - loc.bit_offset));
                }
            }
        }
    }
}

// specialization to pack from 32 bits to basetype of 32 bits
template<>
void pack<uint32_t, uint32_t>(uint32_t* buf, uint32_t* packed_data, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("32bit packing\n");
    const int bits_type = 32;
    int col;
    // main loop
    for (col = 0; col < col_num - 127; col+=128) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1 = _mm512_loadu_si512((void*)(packed_loc_init + 16));
        auto p_vec_2 = _mm512_loadu_si512((void*)(packed_loc_init + 32));
        auto p_vec_3 = _mm512_loadu_si512((void*)(packed_loc_init + 48));
        auto p_vec_4 = _mm512_loadu_si512((void*)(packed_loc_init + 64));
        auto p_vec_5 = _mm512_loadu_si512((void*)(packed_loc_init + 80));
        auto p_vec_6 = _mm512_loadu_si512((void*)(packed_loc_init + 96));
        auto p_vec_7 = _mm512_loadu_si512((void*)(packed_loc_init + 112));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint32_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));
            auto d_vec_1 = _mm512_loadu_si512((void*)(buf_loc + 16));
            auto d_vec_2 = _mm512_loadu_si512((void*)(buf_loc + 32));
            auto d_vec_3 = _mm512_loadu_si512((void*)(buf_loc + 48));
            auto d_vec_4 = _mm512_loadu_si512((void*)(buf_loc + 64));
            auto d_vec_5 = _mm512_loadu_si512((void*)(buf_loc + 80));
            auto d_vec_6 = _mm512_loadu_si512((void*)(buf_loc + 96));
            auto d_vec_7 = _mm512_loadu_si512((void*)(buf_loc + 112));

            // shift data vectors to right spot
            auto sl_d_vec_0 = _mm512_slli_epi32(d_vec_0, col_loc.bit_offset);
            auto sl_d_vec_1 = _mm512_slli_epi32(d_vec_1, col_loc.bit_offset);
            auto sl_d_vec_2 = _mm512_slli_epi32(d_vec_2, col_loc.bit_offset);
            auto sl_d_vec_3 = _mm512_slli_epi32(d_vec_3, col_loc.bit_offset);
            auto sl_d_vec_4 = _mm512_slli_epi32(d_vec_4, col_loc.bit_offset);
            auto sl_d_vec_5 = _mm512_slli_epi32(d_vec_5, col_loc.bit_offset);
            auto sl_d_vec_6 = _mm512_slli_epi32(d_vec_6, col_loc.bit_offset);
            auto sl_d_vec_7 = _mm512_slli_epi32(d_vec_7, col_loc.bit_offset);

            // or
            p_vec_0 = _mm512_or_si512(p_vec_0, sl_d_vec_0);
            p_vec_1 = _mm512_or_si512(p_vec_1, sl_d_vec_1);
            p_vec_2 = _mm512_or_si512(p_vec_2, sl_d_vec_2);
            p_vec_3 = _mm512_or_si512(p_vec_3, sl_d_vec_3);
            p_vec_4 = _mm512_or_si512(p_vec_4, sl_d_vec_4);
            p_vec_5 = _mm512_or_si512(p_vec_5, sl_d_vec_5);
            p_vec_6 = _mm512_or_si512(p_vec_6, sl_d_vec_6);
            p_vec_7 = _mm512_or_si512(p_vec_7, sl_d_vec_7);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0);
                _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1);
                _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2);
                _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3);
                _mm512_storeu_si512((void*)(store_loc + 64), p_vec_4);
                _mm512_storeu_si512((void*)(store_loc + 80), p_vec_5);
                _mm512_storeu_si512((void*)(store_loc + 96), p_vec_6);
                _mm512_storeu_si512((void*)(store_loc + 112), p_vec_7);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0 = _mm512_srli_epi32(d_vec_0, shift_amount);
                    p_vec_1 = _mm512_srli_epi32(d_vec_1, shift_amount);
                    p_vec_2 = _mm512_srli_epi32(d_vec_2, shift_amount);
                    p_vec_3 = _mm512_srli_epi32(d_vec_3, shift_amount);
                    p_vec_4 = _mm512_srli_epi32(d_vec_4, shift_amount);
                    p_vec_5 = _mm512_srli_epi32(d_vec_5, shift_amount);
                    p_vec_6 = _mm512_srli_epi32(d_vec_6, shift_amount);
                    p_vec_7 = _mm512_srli_epi32(d_vec_7, shift_amount);
                } else {
                    p_vec_0 = _mm512_setzero_si512();
                    p_vec_1 = _mm512_setzero_si512();
                    p_vec_2 = _mm512_setzero_si512();
                    p_vec_3 = _mm512_setzero_si512();
                    p_vec_4 = _mm512_setzero_si512();
                    p_vec_5 = _mm512_setzero_si512();
                    p_vec_6 = _mm512_setzero_si512();
                    p_vec_7 = _mm512_setzero_si512();
                }
            }
        }

        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0);
        _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1);
        _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2);
        _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3);
        _mm512_storeu_si512((void*)(store_loc + 64), p_vec_4);
        _mm512_storeu_si512((void*)(store_loc + 80), p_vec_5);
        _mm512_storeu_si512((void*)(store_loc + 96), p_vec_6);
        _mm512_storeu_si512((void*)(store_loc + 112), p_vec_7);
    }

    // first cleanup
    for (; col < col_num - 15; col+=16) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint32_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));

            // shift data vectors to right spot
            auto sl_d_vec_0 = _mm512_slli_epi32(d_vec_0, col_loc.bit_offset);

            // or
            p_vec_0 = _mm512_or_si512(p_vec_0, sl_d_vec_0);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0 = _mm512_srli_epi32(d_vec_0, shift_amount);
                } else {
                    p_vec_0 = _mm512_setzero_si512();
                }
            }
        }
        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0);
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        // go through columns
        for (int col1 = col; col1 < col_num; col1++) {
            // store bits
            packed_data[loc.packed_num + col1] |= (uint32_t)(buf[row * buf_cols + col1] << loc.bit_offset);
        }

        // adjust current bit position
        loc.bit_offset += num_bits;

        // wrap around
        if (loc.bit_offset >= bits_type) {
            loc.bit_offset &= (bits_type - 1);
            loc.packed_num += loc.eff_cols;

            // initialize for next run
            if (loc.bit_offset != 0) {
                for (int col1 = col; col1 < col_num; col1++) {
                    // store remaining bits from last element
                    packed_data[loc.packed_num + col1] = (uint32_t)(buf[row * buf_cols + col1] >> (num_bits - loc.bit_offset));
                }
            }
        }
    }
}

// specialization to pack from 16 bits to basetype of 32 bits
template<>
void pack<uint16_t, uint32_t>(uint16_t* buf, uint32_t* packed_data, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("16bit packing\n");
    const int bits_type = 32;
    int col;
    auto permute_vec = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    // main loop
    for (col = 0; col < col_num - 127; col+=128) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc_init + 16));

        auto p_vec_0_1 = _mm512_loadu_si512((void*)(packed_loc_init + 32));
        auto p_vec_1_1 = _mm512_loadu_si512((void*)(packed_loc_init + 48));

        auto p_vec_0_2 = _mm512_loadu_si512((void*)(packed_loc_init + 64));
        auto p_vec_1_2 = _mm512_loadu_si512((void*)(packed_loc_init + 80));

        auto p_vec_0_3 = _mm512_loadu_si512((void*)(packed_loc_init + 96));
        auto p_vec_1_3 = _mm512_loadu_si512((void*)(packed_loc_init + 112));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint16_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));
            auto d_vec_left_0 = _mm512_castsi512_si256(d_vec_0);
            auto d_vec_right_0 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_0));

            auto d_vec_1 = _mm512_loadu_si512((void*)(buf_loc + 32));
            auto d_vec_left_1 = _mm512_castsi512_si256(d_vec_1);
            auto d_vec_right_1 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_1));

            auto d_vec_2 = _mm512_loadu_si512((void*)(buf_loc + 64));
            auto d_vec_left_2 = _mm512_castsi512_si256(d_vec_2);
            auto d_vec_right_2 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_2));

            auto d_vec_3 = _mm512_loadu_si512((void*)(buf_loc + 96));
            auto d_vec_left_3 = _mm512_castsi512_si256(d_vec_3);
            auto d_vec_right_3 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_3));

            auto scaled_d_vec_0_0 = _mm512_cvtepu16_epi32(d_vec_left_0);
            auto scaled_d_vec_1_0 = _mm512_cvtepu16_epi32(d_vec_right_0);

            auto scaled_d_vec_0_1 = _mm512_cvtepu16_epi32(d_vec_left_1);
            auto scaled_d_vec_1_1 = _mm512_cvtepu16_epi32(d_vec_right_1);

            auto scaled_d_vec_0_2 = _mm512_cvtepu16_epi32(d_vec_left_2);
            auto scaled_d_vec_1_2 = _mm512_cvtepu16_epi32(d_vec_right_2);

            auto scaled_d_vec_0_3 = _mm512_cvtepu16_epi32(d_vec_left_3);
            auto scaled_d_vec_1_3 = _mm512_cvtepu16_epi32(d_vec_right_3);

            // shift data vectors to right spot
            auto sl_d_vec_0_0 = _mm512_slli_epi32(scaled_d_vec_0_0, col_loc.bit_offset);
            auto sl_d_vec_1_0 = _mm512_slli_epi32(scaled_d_vec_1_0, col_loc.bit_offset);

            auto sl_d_vec_0_1 = _mm512_slli_epi32(scaled_d_vec_0_1, col_loc.bit_offset);
            auto sl_d_vec_1_1 = _mm512_slli_epi32(scaled_d_vec_1_1, col_loc.bit_offset);

            auto sl_d_vec_0_2 = _mm512_slli_epi32(scaled_d_vec_0_2, col_loc.bit_offset);
            auto sl_d_vec_1_2 = _mm512_slli_epi32(scaled_d_vec_1_2, col_loc.bit_offset);

            auto sl_d_vec_0_3 = _mm512_slli_epi32(scaled_d_vec_0_3, col_loc.bit_offset);
            auto sl_d_vec_1_3 = _mm512_slli_epi32(scaled_d_vec_1_3, col_loc.bit_offset);

            // or
            p_vec_0_0 = _mm512_or_si512(p_vec_0_0, sl_d_vec_0_0);
            p_vec_1_0 = _mm512_or_si512(p_vec_1_0, sl_d_vec_1_0);

            p_vec_0_1 = _mm512_or_si512(p_vec_0_1, sl_d_vec_0_1);
            p_vec_1_1 = _mm512_or_si512(p_vec_1_1, sl_d_vec_1_1);

            p_vec_0_2 = _mm512_or_si512(p_vec_0_2, sl_d_vec_0_2);
            p_vec_1_2 = _mm512_or_si512(p_vec_1_2, sl_d_vec_1_2);

            p_vec_0_3 = _mm512_or_si512(p_vec_0_3, sl_d_vec_0_3);
            p_vec_1_3 = _mm512_or_si512(p_vec_1_3, sl_d_vec_1_3);
            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
                _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);

                _mm512_storeu_si512((void*)(store_loc + 32), p_vec_0_1);
                _mm512_storeu_si512((void*)(store_loc + 48), p_vec_1_1);

                _mm512_storeu_si512((void*)(store_loc + 64), p_vec_0_2);
                _mm512_storeu_si512((void*)(store_loc + 80), p_vec_1_2);

                _mm512_storeu_si512((void*)(store_loc + 96), p_vec_0_3);
                _mm512_storeu_si512((void*)(store_loc + 112), p_vec_1_3);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0_0 = _mm512_srli_epi32(scaled_d_vec_0_0, shift_amount);
                    p_vec_1_0 = _mm512_srli_epi32(scaled_d_vec_1_0, shift_amount);

                    p_vec_0_1 = _mm512_srli_epi32(scaled_d_vec_0_1, shift_amount);
                    p_vec_1_1 = _mm512_srli_epi32(scaled_d_vec_1_1, shift_amount);

                    p_vec_0_2 = _mm512_srli_epi32(scaled_d_vec_0_2, shift_amount);
                    p_vec_1_2 = _mm512_srli_epi32(scaled_d_vec_1_2, shift_amount);

                    p_vec_0_3 = _mm512_srli_epi32(scaled_d_vec_0_3, shift_amount);
                    p_vec_1_3 = _mm512_srli_epi32(scaled_d_vec_1_3, shift_amount);
                } else {
                    p_vec_0_0 = _mm512_setzero_si512();
                    p_vec_1_0 = _mm512_setzero_si512();

                    p_vec_0_1 = _mm512_setzero_si512();
                    p_vec_1_1 = _mm512_setzero_si512();

                    p_vec_0_2 = _mm512_setzero_si512();
                    p_vec_1_2 = _mm512_setzero_si512();

                    p_vec_0_3 = _mm512_setzero_si512();
                    p_vec_1_3 = _mm512_setzero_si512();
                }
            }
        }

        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
        _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);

        _mm512_storeu_si512((void*)(store_loc + 32), p_vec_0_1);
        _mm512_storeu_si512((void*)(store_loc + 48), p_vec_1_1);

        _mm512_storeu_si512((void*)(store_loc + 64), p_vec_0_2);
        _mm512_storeu_si512((void*)(store_loc + 80), p_vec_1_2);

        _mm512_storeu_si512((void*)(store_loc + 96), p_vec_0_3);
        _mm512_storeu_si512((void*)(store_loc + 112), p_vec_1_3);
    }

    // cleanup
    for (; col < col_num - 31; col+=32) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc_init + 16));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint16_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));
            auto d_vec_left_0 = _mm512_castsi512_si256(d_vec_0);
            auto d_vec_right_0 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_0));

            auto scaled_d_vec_0_0 = _mm512_cvtepu16_epi32(d_vec_left_0);
            auto scaled_d_vec_1_0 = _mm512_cvtepu16_epi32(d_vec_right_0);

            // shift data vectors to right spot
            auto sl_d_vec_0_0 = _mm512_slli_epi32(scaled_d_vec_0_0, col_loc.bit_offset);
            auto sl_d_vec_1_0 = _mm512_slli_epi32(scaled_d_vec_1_0, col_loc.bit_offset);

            // or
            p_vec_0_0 = _mm512_or_si512(p_vec_0_0, sl_d_vec_0_0);
            p_vec_1_0 = _mm512_or_si512(p_vec_1_0, sl_d_vec_1_0);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
                _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0_0 = _mm512_srli_epi32(scaled_d_vec_0_0, shift_amount);
                    p_vec_1_0 = _mm512_srli_epi32(scaled_d_vec_1_0, shift_amount);
                } else {
                    p_vec_0_0 = _mm512_setzero_si512();
                    p_vec_1_0 = _mm512_setzero_si512();
                }
            }
        }
        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
        _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        // go through columns
        for (int col1 = col; col1 < col_num; col1++) {
            // store bits
            packed_data[loc.packed_num + col1] |= (uint32_t)(buf[row * buf_cols + col1] << loc.bit_offset);
        }

        // adjust current bit position
        loc.bit_offset += num_bits;

        // wrap around
        if (loc.bit_offset >= bits_type) {
            loc.bit_offset &= (bits_type - 1);
            loc.packed_num += loc.eff_cols;

            // initialize for next run
            if (loc.bit_offset != 0) {
                for (int col1 = col; col1 < col_num; col1++) {
                    // store remaining bits from last element
                    packed_data[loc.packed_num + col1] = (uint32_t)(buf[row * buf_cols + col1] >> (num_bits - loc.bit_offset));
                }
            }
        }
    }
}

// specialization to pack from 8 bits to basetype of 32 bits
template<>
void pack<uint8_t, uint32_t>(uint8_t* buf, uint32_t* packed_data, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("8bit packing\n");
    const int bits_type = 32;
    int col;
    auto permute_vec_1 = _mm512_set_epi64(7, 6, 5, 4, 1, 0, 3, 2);
    auto permute_vec_2 = _mm512_set_epi64(7, 6, 1, 0, 3, 2, 5, 4);
    auto permute_vec_3 = _mm512_set_epi64(1, 0, 5, 4, 3, 2, 7, 6);
    // main loop
    for (col = 0; col < col_num - 127; col+=128) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc_init + 16));
        auto p_vec_2_0 = _mm512_loadu_si512((void*)(packed_loc_init + 32));
        auto p_vec_3_0 = _mm512_loadu_si512((void*)(packed_loc_init + 48));

        auto p_vec_0_1 = _mm512_loadu_si512((void*)(packed_loc_init + 64));
        auto p_vec_1_1 = _mm512_loadu_si512((void*)(packed_loc_init + 80));
        auto p_vec_2_1 = _mm512_loadu_si512((void*)(packed_loc_init + 96));
        auto p_vec_3_1 = _mm512_loadu_si512((void*)(packed_loc_init + 112));
        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint8_t* buf_loc = buf + row*buf_cols + col;
            // |0|1|2|3|
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));
            auto d_vec_1 = _mm512_loadu_si512((void*)(buf_loc + 64));
            // 0|...
            auto d_vec_0_0 = _mm512_castsi512_si128(d_vec_0);
            auto d_vec_0_1 = _mm512_castsi512_si128(d_vec_1);
            // 1|...
            auto d_vec_1_0 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_1, d_vec_0));
            auto d_vec_1_1 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_1, d_vec_1));
            // 2|...
            auto d_vec_2_0 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_2, d_vec_0));
            auto d_vec_2_1 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_2, d_vec_1));
            // 3|...
            auto d_vec_3_0 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_3, d_vec_0));
            auto d_vec_3_1 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_3, d_vec_1));

            auto scaled_d_vec_0_0 = _mm512_cvtepu8_epi32(d_vec_0_0);
            auto scaled_d_vec_1_0 = _mm512_cvtepu8_epi32(d_vec_1_0);
            auto scaled_d_vec_2_0 = _mm512_cvtepu8_epi32(d_vec_2_0);
            auto scaled_d_vec_3_0 = _mm512_cvtepu8_epi32(d_vec_3_0);

            auto scaled_d_vec_0_1 = _mm512_cvtepu8_epi32(d_vec_0_1);
            auto scaled_d_vec_1_1 = _mm512_cvtepu8_epi32(d_vec_1_1);
            auto scaled_d_vec_2_1 = _mm512_cvtepu8_epi32(d_vec_2_1);
            auto scaled_d_vec_3_1 = _mm512_cvtepu8_epi32(d_vec_3_1);

            // shift data vectors to right spot
            auto sl_d_vec_0_0 = _mm512_slli_epi32(scaled_d_vec_0_0, col_loc.bit_offset);
            auto sl_d_vec_1_0 = _mm512_slli_epi32(scaled_d_vec_1_0, col_loc.bit_offset);
            auto sl_d_vec_2_0 = _mm512_slli_epi32(scaled_d_vec_2_0, col_loc.bit_offset);
            auto sl_d_vec_3_0 = _mm512_slli_epi32(scaled_d_vec_3_0, col_loc.bit_offset);

            auto sl_d_vec_0_1 = _mm512_slli_epi32(scaled_d_vec_0_1, col_loc.bit_offset);
            auto sl_d_vec_1_1 = _mm512_slli_epi32(scaled_d_vec_1_1, col_loc.bit_offset);
            auto sl_d_vec_2_1 = _mm512_slli_epi32(scaled_d_vec_2_1, col_loc.bit_offset);
            auto sl_d_vec_3_1 = _mm512_slli_epi32(scaled_d_vec_3_1, col_loc.bit_offset);

            // or
            p_vec_0_0 = _mm512_or_si512(p_vec_0_0, sl_d_vec_0_0);
            p_vec_1_0 = _mm512_or_si512(p_vec_1_0, sl_d_vec_1_0);
            p_vec_2_0 = _mm512_or_si512(p_vec_2_0, sl_d_vec_2_0);
            p_vec_3_0 = _mm512_or_si512(p_vec_3_0, sl_d_vec_3_0);

            p_vec_0_1 = _mm512_or_si512(p_vec_0_1, sl_d_vec_0_1);
            p_vec_1_1 = _mm512_or_si512(p_vec_1_1, sl_d_vec_1_1);
            p_vec_2_1 = _mm512_or_si512(p_vec_2_1, sl_d_vec_2_1);
            p_vec_3_1 = _mm512_or_si512(p_vec_3_1, sl_d_vec_3_1);
            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
                _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);
                _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2_0);
                _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3_0);

                _mm512_storeu_si512((void*)(store_loc + 64), p_vec_0_1);
                _mm512_storeu_si512((void*)(store_loc + 80), p_vec_1_1);
                _mm512_storeu_si512((void*)(store_loc + 96), p_vec_2_1);
                _mm512_storeu_si512((void*)(store_loc + 112), p_vec_3_1);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0_0 = _mm512_srli_epi32(scaled_d_vec_0_0, shift_amount);
                    p_vec_1_0 = _mm512_srli_epi32(scaled_d_vec_1_0, shift_amount);
                    p_vec_2_0 = _mm512_srli_epi32(scaled_d_vec_2_0, shift_amount);
                    p_vec_3_0 = _mm512_srli_epi32(scaled_d_vec_3_0, shift_amount);

                    p_vec_0_1 = _mm512_srli_epi32(scaled_d_vec_0_1, shift_amount);
                    p_vec_1_1 = _mm512_srli_epi32(scaled_d_vec_1_1, shift_amount);
                    p_vec_2_1 = _mm512_srli_epi32(scaled_d_vec_2_1, shift_amount);
                    p_vec_3_1 = _mm512_srli_epi32(scaled_d_vec_3_1, shift_amount);
                } else {
                    p_vec_0_0 = _mm512_setzero_si512();
                    p_vec_1_0 = _mm512_setzero_si512();
                    p_vec_2_0 = _mm512_setzero_si512();
                    p_vec_3_0 = _mm512_setzero_si512();

                    p_vec_0_1 = _mm512_setzero_si512();
                    p_vec_1_1 = _mm512_setzero_si512();
                    p_vec_2_1 = _mm512_setzero_si512();
                    p_vec_3_1 = _mm512_setzero_si512();
                }
            }
        }
        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
        _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);
        _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2_0);
        _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3_0);

        _mm512_storeu_si512((void*)(store_loc + 64), p_vec_0_1);
        _mm512_storeu_si512((void*)(store_loc + 80), p_vec_1_1);
        _mm512_storeu_si512((void*)(store_loc + 96), p_vec_2_1);
        _mm512_storeu_si512((void*)(store_loc + 112), p_vec_3_1);
    }

    // cleanup
    for (; col < col_num - 63; col+=64) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc_init + 16));
        auto p_vec_2_0 = _mm512_loadu_si512((void*)(packed_loc_init + 32));
        auto p_vec_3_0 = _mm512_loadu_si512((void*)(packed_loc_init + 48));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint8_t* buf_loc = buf + row*buf_cols + col;
            // |0|1|2|3|
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));
            // 0|...
            auto d_vec_0_0 = _mm512_castsi512_si128(d_vec_0);
            // 1|...
            auto d_vec_1_0 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_1, d_vec_0));
            // 2|...
            auto d_vec_2_0 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_2, d_vec_0));
            // 3|...
            auto d_vec_3_0 = _mm512_castsi512_si128(_mm512_permutexvar_epi64(permute_vec_3, d_vec_0));

            auto scaled_d_vec_0_0 = _mm512_cvtepu8_epi32(d_vec_0_0);
            auto scaled_d_vec_1_0 = _mm512_cvtepu8_epi32(d_vec_1_0);
            auto scaled_d_vec_2_0 = _mm512_cvtepu8_epi32(d_vec_2_0);
            auto scaled_d_vec_3_0 = _mm512_cvtepu8_epi32(d_vec_3_0);

            // shift data vectors to right spot
            auto sl_d_vec_0_0 = _mm512_slli_epi32(scaled_d_vec_0_0, col_loc.bit_offset);
            auto sl_d_vec_1_0 = _mm512_slli_epi32(scaled_d_vec_1_0, col_loc.bit_offset);
            auto sl_d_vec_2_0 = _mm512_slli_epi32(scaled_d_vec_2_0, col_loc.bit_offset);
            auto sl_d_vec_3_0 = _mm512_slli_epi32(scaled_d_vec_3_0, col_loc.bit_offset);

            // or
            p_vec_0_0 = _mm512_or_si512(p_vec_0_0, sl_d_vec_0_0);
            p_vec_1_0 = _mm512_or_si512(p_vec_1_0, sl_d_vec_1_0);
            p_vec_2_0 = _mm512_or_si512(p_vec_2_0, sl_d_vec_2_0);
            p_vec_3_0 = _mm512_or_si512(p_vec_3_0, sl_d_vec_3_0);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
                _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);
                _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2_0);
                _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3_0);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0_0 = _mm512_srli_epi32(scaled_d_vec_0_0, shift_amount);
                    p_vec_1_0 = _mm512_srli_epi32(scaled_d_vec_1_0, shift_amount);
                    p_vec_2_0 = _mm512_srli_epi32(scaled_d_vec_2_0, shift_amount);
                    p_vec_3_0 = _mm512_srli_epi32(scaled_d_vec_3_0, shift_amount);
                } else {
                    p_vec_0_0 = _mm512_setzero_si512();
                    p_vec_1_0 = _mm512_setzero_si512();
                    p_vec_2_0 = _mm512_setzero_si512();
                    p_vec_3_0 = _mm512_setzero_si512();
                }
            }
        }

        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0_0);
        _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1_0);
        _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2_0);
        _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3_0);
    }
    // final cleanup
    for (int row = 0; row < row_num; row++) {
        // go through columns
        for (int col1 = col; col1 < col_num; col1++) {
            // store bits
            packed_data[loc.packed_num + col1] |= (uint32_t)(buf[row * buf_cols + col1] << loc.bit_offset);
        }

        // adjust current bit position
        loc.bit_offset += num_bits;

        // wrap around
        if (loc.bit_offset >= bits_type) {
            loc.bit_offset &= (bits_type -1);
            loc.packed_num += loc.eff_cols;

            // initialize for next run
            if (loc.bit_offset != 0) {
                for (int col1 = col; col1 < col_num; col1++) {
                    // store remaining bits from last element
                    packed_data[loc.packed_num + col1] = (uint32_t)(buf[row * buf_cols + col1] >> (num_bits - loc.bit_offset));
                }
            }
        }
    }
}

// specialization to pack from 64 bits to basetype of 32 bits
template<>
void pack<uint64_t, uint32_t>(uint64_t* buf, uint32_t* packed_data, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("64bit packing\n");
    const int bits_type = 32;
    int col;
    // main loop
    for (col = 0; col < col_num - 63; col+=64) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1 = _mm512_loadu_si512((void*)(packed_loc_init + 16));
        auto p_vec_2 = _mm512_loadu_si512((void*)(packed_loc_init + 32));
        auto p_vec_3 = _mm512_loadu_si512((void*)(packed_loc_init + 48));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint64_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0_0 = _mm512_loadu_si512((void*)(buf_loc));
            auto d_vec_1_0 = _mm512_loadu_si512((void*)(buf_loc + 8));
            auto d_vec_0_1 = _mm512_loadu_si512((void*)(buf_loc + 16));
            auto d_vec_1_1 = _mm512_loadu_si512((void*)(buf_loc + 24));
            auto d_vec_0_2 = _mm512_loadu_si512((void*)(buf_loc + 32));
            auto d_vec_1_2 = _mm512_loadu_si512((void*)(buf_loc + 40));
            auto d_vec_0_3 = _mm512_loadu_si512((void*)(buf_loc + 48));
            auto d_vec_1_3 = _mm512_loadu_si512((void*)(buf_loc + 56));

            // convert to 32 bit ints
            auto cvt_d_vec_0_0 = _mm512_cvtepi64_epi32(d_vec_0_0);
            auto cvt_d_vec_1_0 = _mm512_cvtepi64_epi32(d_vec_1_0);
            auto cvt_d_vec_0_1 = _mm512_cvtepi64_epi32(d_vec_0_1);
            auto cvt_d_vec_1_1 = _mm512_cvtepi64_epi32(d_vec_1_1);
            auto cvt_d_vec_0_2 = _mm512_cvtepi64_epi32(d_vec_0_2);
            auto cvt_d_vec_1_2 = _mm512_cvtepi64_epi32(d_vec_1_2);
            auto cvt_d_vec_0_3 = _mm512_cvtepi64_epi32(d_vec_0_3);
            auto cvt_d_vec_1_3 = _mm512_cvtepi64_epi32(d_vec_1_3);

            // put into one register
            auto d_vec_0 = _mm512_inserti64x4(_mm512_castsi256_si512(cvt_d_vec_0_0), cvt_d_vec_1_0, 1);
            auto d_vec_1 = _mm512_inserti64x4(_mm512_castsi256_si512(cvt_d_vec_0_1), cvt_d_vec_1_1, 1);
            auto d_vec_2 = _mm512_inserti64x4(_mm512_castsi256_si512(cvt_d_vec_0_2), cvt_d_vec_1_2, 1);
            auto d_vec_3 = _mm512_inserti64x4(_mm512_castsi256_si512(cvt_d_vec_0_3), cvt_d_vec_1_3, 1);

            // shift data vectors to right spot
            auto sl_d_vec_0 = _mm512_slli_epi32(d_vec_0, col_loc.bit_offset);
            auto sl_d_vec_1 = _mm512_slli_epi32(d_vec_1, col_loc.bit_offset);
            auto sl_d_vec_2 = _mm512_slli_epi32(d_vec_2, col_loc.bit_offset);
            auto sl_d_vec_3 = _mm512_slli_epi32(d_vec_3, col_loc.bit_offset);

            // or
            p_vec_0 = _mm512_or_si512(p_vec_0, sl_d_vec_0);
            p_vec_1 = _mm512_or_si512(p_vec_1, sl_d_vec_1);
            p_vec_2 = _mm512_or_si512(p_vec_2, sl_d_vec_2);
            p_vec_3 = _mm512_or_si512(p_vec_3, sl_d_vec_3);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0);
                _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1);
                _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2);
                _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0 = _mm512_srli_epi32(d_vec_0, shift_amount);
                    p_vec_1 = _mm512_srli_epi32(d_vec_1, shift_amount);
                    p_vec_2 = _mm512_srli_epi32(d_vec_2, shift_amount);
                    p_vec_3 = _mm512_srli_epi32(d_vec_3, shift_amount);
                } else {
                    p_vec_0 = _mm512_setzero_si512();
                    p_vec_1 = _mm512_setzero_si512();
                    p_vec_2 = _mm512_setzero_si512();
                    p_vec_3 = _mm512_setzero_si512();
                }
            }
        }

        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0);
        _mm512_storeu_si512((void*)(store_loc + 16), p_vec_1);
        _mm512_storeu_si512((void*)(store_loc + 32), p_vec_2);
        _mm512_storeu_si512((void*)(store_loc + 48), p_vec_3);
    }

    // first cleanup
    for (; col < col_num - 15; col+=16) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint64_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0_0 = _mm512_loadu_si512((void*)(buf_loc));
            auto d_vec_1_0 = _mm512_loadu_si512((void*)(buf_loc + 8));

            auto cvt_d_vec_0_0 = _mm512_cvtepi64_epi32(d_vec_0_0);
            auto cvt_d_vec_1_0 = _mm512_cvtepi64_epi32(d_vec_1_0);

            auto d_vec_0 = _mm512_inserti64x4(_mm512_castsi256_si512(cvt_d_vec_0_0), cvt_d_vec_1_0, 1);

            // shift data vectors to right spot
            auto sl_d_vec_0 = _mm512_slli_epi32(d_vec_0, col_loc.bit_offset);

            // or
            p_vec_0 = _mm512_or_si512(p_vec_0, sl_d_vec_0);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint32_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0 = _mm512_srli_epi32(d_vec_0, shift_amount);
                } else {
                    p_vec_0 = _mm512_setzero_si512();
                }
            }
        }
        // store
        uint32_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0);
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        // go through columns
        for (int col1 = col; col1 < col_num; col1++) {
            // store bits
            packed_data[loc.packed_num + col1] |= (uint32_t)(buf[row * buf_cols + col1] << loc.bit_offset);
        }

        // adjust current bit position
        loc.bit_offset += num_bits;

        // wrap around
        if (loc.bit_offset >= bits_type) {
            loc.bit_offset &= (bits_type - 1);
            loc.packed_num += loc.eff_cols;

            // initialize for next run
            if (loc.bit_offset != 0) {
                for (int col1 = col; col1 < col_num; col1++) {
                    // store remaining bits from last element
                    packed_data[loc.packed_num + col1] = (uint32_t)(buf[row * buf_cols + col1] >> (num_bits - loc.bit_offset));
                }
            }
        }
    }
}

// specialization to pack from 64 bits to basetype of 64 bits
template<>
void pack<uint64_t, uint64_t>(uint64_t* buf, uint64_t* packed_data, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("64bit to 64bit packing\n");
    const int bits_type = 64;
    int col;
    // main loop
    for (col = 0; col < col_num - 63; col+=64) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint64_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1 = _mm512_loadu_si512((void*)(packed_loc_init + 8));
        auto p_vec_2 = _mm512_loadu_si512((void*)(packed_loc_init + 16));
        auto p_vec_3 = _mm512_loadu_si512((void*)(packed_loc_init + 24));
        auto p_vec_4 = _mm512_loadu_si512((void*)(packed_loc_init + 32));
        auto p_vec_5 = _mm512_loadu_si512((void*)(packed_loc_init + 40));
        auto p_vec_6 = _mm512_loadu_si512((void*)(packed_loc_init + 48));
        auto p_vec_7 = _mm512_loadu_si512((void*)(packed_loc_init + 56));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint64_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));
            auto d_vec_1 = _mm512_loadu_si512((void*)(buf_loc + 8));
            auto d_vec_2 = _mm512_loadu_si512((void*)(buf_loc + 16));
            auto d_vec_3 = _mm512_loadu_si512((void*)(buf_loc + 24));
            auto d_vec_4 = _mm512_loadu_si512((void*)(buf_loc + 32));
            auto d_vec_5 = _mm512_loadu_si512((void*)(buf_loc + 40));
            auto d_vec_6 = _mm512_loadu_si512((void*)(buf_loc + 48));
            auto d_vec_7 = _mm512_loadu_si512((void*)(buf_loc + 56));

            // shift data vectors to right spot
            auto sl_d_vec_0 = _mm512_slli_epi64(d_vec_0, col_loc.bit_offset);
            auto sl_d_vec_1 = _mm512_slli_epi64(d_vec_1, col_loc.bit_offset);
            auto sl_d_vec_2 = _mm512_slli_epi64(d_vec_2, col_loc.bit_offset);
            auto sl_d_vec_3 = _mm512_slli_epi64(d_vec_3, col_loc.bit_offset);
            auto sl_d_vec_4 = _mm512_slli_epi64(d_vec_4, col_loc.bit_offset);
            auto sl_d_vec_5 = _mm512_slli_epi64(d_vec_5, col_loc.bit_offset);
            auto sl_d_vec_6 = _mm512_slli_epi64(d_vec_6, col_loc.bit_offset);
            auto sl_d_vec_7 = _mm512_slli_epi64(d_vec_7, col_loc.bit_offset);

            // or
            p_vec_0 = _mm512_or_si512(p_vec_0, sl_d_vec_0);
            p_vec_1 = _mm512_or_si512(p_vec_1, sl_d_vec_1);
            p_vec_2 = _mm512_or_si512(p_vec_2, sl_d_vec_2);
            p_vec_3 = _mm512_or_si512(p_vec_3, sl_d_vec_3);
            p_vec_4 = _mm512_or_si512(p_vec_4, sl_d_vec_4);
            p_vec_5 = _mm512_or_si512(p_vec_5, sl_d_vec_5);
            p_vec_6 = _mm512_or_si512(p_vec_6, sl_d_vec_6);
            p_vec_7 = _mm512_or_si512(p_vec_7, sl_d_vec_7);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint64_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0);
                _mm512_storeu_si512((void*)(store_loc + 8), p_vec_1);
                _mm512_storeu_si512((void*)(store_loc + 16), p_vec_2);
                _mm512_storeu_si512((void*)(store_loc + 24), p_vec_3);
                _mm512_storeu_si512((void*)(store_loc + 32), p_vec_4);
                _mm512_storeu_si512((void*)(store_loc + 40), p_vec_5);
                _mm512_storeu_si512((void*)(store_loc + 48), p_vec_6);
                _mm512_storeu_si512((void*)(store_loc + 56), p_vec_7);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0 = _mm512_srli_epi64(d_vec_0, shift_amount);
                    p_vec_1 = _mm512_srli_epi64(d_vec_1, shift_amount);
                    p_vec_2 = _mm512_srli_epi64(d_vec_2, shift_amount);
                    p_vec_3 = _mm512_srli_epi64(d_vec_3, shift_amount);
                    p_vec_4 = _mm512_srli_epi64(d_vec_4, shift_amount);
                    p_vec_5 = _mm512_srli_epi64(d_vec_5, shift_amount);
                    p_vec_6 = _mm512_srli_epi64(d_vec_6, shift_amount);
                    p_vec_7 = _mm512_srli_epi64(d_vec_7, shift_amount);
                } else {
                    p_vec_0 = _mm512_setzero_si512();
                    p_vec_1 = _mm512_setzero_si512();
                    p_vec_2 = _mm512_setzero_si512();
                    p_vec_3 = _mm512_setzero_si512();
                    p_vec_4 = _mm512_setzero_si512();
                    p_vec_5 = _mm512_setzero_si512();
                    p_vec_6 = _mm512_setzero_si512();
                    p_vec_7 = _mm512_setzero_si512();
                }
            }
        }

        // store
        uint64_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0);
        _mm512_storeu_si512((void*)(store_loc + 8), p_vec_1);
        _mm512_storeu_si512((void*)(store_loc + 16), p_vec_2);
        _mm512_storeu_si512((void*)(store_loc + 24), p_vec_3);
        _mm512_storeu_si512((void*)(store_loc + 32), p_vec_4);
        _mm512_storeu_si512((void*)(store_loc + 40), p_vec_5);
        _mm512_storeu_si512((void*)(store_loc + 48), p_vec_6);
        _mm512_storeu_si512((void*)(store_loc + 56), p_vec_7);
    }

    // first cleanup
    for (; col < col_num - 7; col+=8) {
        Packed_Loc col_loc = loc;
        // load packed vectors
        uint64_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        for (int row = 0; row < row_num; row++) {
            // load data vectors
            uint64_t* buf_loc = buf + row*buf_cols + col;
            auto d_vec_0 = _mm512_loadu_si512((void*)(buf_loc));

            // shift data vectors to right spot
            auto sl_d_vec_0 = _mm512_slli_epi64(d_vec_0, col_loc.bit_offset);

            // or
            p_vec_0 = _mm512_or_si512(p_vec_0, sl_d_vec_0);

            // adjust current bit position
            col_loc.bit_offset += num_bits;

            // wrap around
            if (col_loc.bit_offset >= bits_type) {
                // store
                uint64_t* store_loc = packed_data + col_loc.packed_num + col;
                _mm512_storeu_si512((void*)(store_loc), p_vec_0);

                col_loc.bit_offset &= (bits_type - 1);
                col_loc.packed_num += col_loc.eff_cols;

                // initialize for next run
                if (col_loc.bit_offset != 0) {
                    int shift_amount = num_bits - col_loc.bit_offset;
                    // store remaining bits from last element
                    p_vec_0 = _mm512_srli_epi64(d_vec_0, shift_amount);
                } else {
                    p_vec_0 = _mm512_setzero_si512();
                }
            }
        }
        // store
        uint64_t* store_loc = packed_data + col_loc.packed_num + col;
        _mm512_storeu_si512((void*)(store_loc), p_vec_0);
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        // go through columns
        for (int col1 = col; col1 < col_num; col1++) {
            // store bits
            packed_data[loc.packed_num + col1] |= (uint64_t)(buf[row * buf_cols + col1] << loc.bit_offset);
        }

        // adjust current bit position
        loc.bit_offset += num_bits;

        // wrap around
        if (loc.bit_offset >= bits_type) {
            loc.bit_offset &= (bits_type - 1);
            loc.packed_num += loc.eff_cols;

            // initialize for next run
            if (loc.bit_offset != 0) {
                for (int col1 = col; col1 < col_num; col1++) {
                    // store remaining bits from last element
                    packed_data[loc.packed_num + col1] = (uint64_t)(buf[row * buf_cols + col1] >> (num_bits - loc.bit_offset));
                }
            }
        }
    }
}
