/**
 * @file APIB_scaled_unpacking.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains implementations of specializations of scaling unpacking functions
 * 
 */

#include "APIB/APIB_scaled_unpacking.h"
#include "APIB/APIB_helper.h"
#include <cstdint>
#include <cstring>

/**
 * @brief Unpacks data in packed format into buffer in format for madd (vnni) and scales each element by alpha
 * 
 * @tparam buf_type Type of destination buffer
 * @tparam base_type Basetype of APIB Mat object containing packed data
 * @param buf Pointer to destination buffer
 * @param packed_data Pointer to packed data
 * @param alpha 
 * @param loc Location of offset
 * @param row_num Number of rows to unpack
 * @param col_num Number of columns to unpack, assumes stays in same block
 * @param buf_cols Number of columns of destination buffer
 * @param num_bits Size of elements in bits
 */
template<typename buf_type, typename base_type>
void scaled_madd_unpack(buf_type* buf, base_type* packed_data, buf_type alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("matrix 32bit unpacking\n");
    uint32_t mask = (uint32_t)((1L << num_bits) - 1); 
    auto mask_vec = _mm512_set1_epi32(mask);
    auto alpha_vec = _mm512_set1_epi16(alpha);
    const int bits_type = 8 * sizeof(base_type);

    int col; 
    // main loop
    for (col=0; col < col_num - 63; col+=64) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        base_type* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1 = _mm512_loadu_si512((void*)(packed_loc_init + 16));
        auto p_vec_2 = _mm512_loadu_si512((void*)(packed_loc_init + 32));
        auto p_vec_3 = _mm512_loadu_si512((void*)(packed_loc_init + 48));

        int row;
        for (row = 0; row < row_num - 1; row+=2) {
            // vec 0
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1_0 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);
            auto sr_p_vec_2_0 = _mm512_srli_epi32(p_vec_2, col_loc.bit_offset);
            auto sr_p_vec_3_0 = _mm512_srli_epi32(p_vec_3, col_loc.bit_offset);

            auto unsc_d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);
            auto unsc_d_vec_1_0 = _mm512_and_si512(mask_vec, sr_p_vec_1_0);
            auto unsc_d_vec_2_0 = _mm512_and_si512(mask_vec, sr_p_vec_2_0);
            auto unsc_d_vec_3_0 = _mm512_and_si512(mask_vec, sr_p_vec_3_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_2 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_3 = _mm512_loadu_si512((void*)(packed_loc + 48));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);
                    auto sl_p_vec_2 = _mm512_slli_epi32(p_vec_2, shift_amount);
                    auto sl_p_vec_3 = _mm512_slli_epi32(p_vec_3, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);
                    auto m_sl_p_vec_2 = _mm512_and_si512(mask_vec, sl_p_vec_2);
                    auto m_sl_p_vec_3 = _mm512_and_si512(mask_vec, sl_p_vec_3);

                    // add to result
                    unsc_d_vec_0_0 = _mm512_or_si512(unsc_d_vec_0_0, m_sl_p_vec_0);
                    unsc_d_vec_1_0 = _mm512_or_si512(unsc_d_vec_1_0, m_sl_p_vec_1);
                    unsc_d_vec_2_0 = _mm512_or_si512(unsc_d_vec_2_0, m_sl_p_vec_2);
                    unsc_d_vec_3_0 = _mm512_or_si512(unsc_d_vec_3_0, m_sl_p_vec_3);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto d_vec_0_0 = _mm512_mullo_epi16(unsc_d_vec_0_0, alpha_vec);
            auto d_vec_1_0 = _mm512_mullo_epi16(unsc_d_vec_1_0, alpha_vec);
            auto d_vec_2_0 = _mm512_mullo_epi16(unsc_d_vec_2_0, alpha_vec);
            auto d_vec_3_0 = _mm512_mullo_epi16(unsc_d_vec_3_0, alpha_vec);

            // vec 1
            // shift and mask
            auto sr_p_vec_0_1 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1_1 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);
            auto sr_p_vec_2_1 = _mm512_srli_epi32(p_vec_2, col_loc.bit_offset);
            auto sr_p_vec_3_1 = _mm512_srli_epi32(p_vec_3, col_loc.bit_offset);

            auto pre_d_vec_0_1 = _mm512_and_si512(mask_vec, sr_p_vec_0_1);
            auto pre_d_vec_1_1 = _mm512_and_si512(mask_vec, sr_p_vec_1_1);
            auto pre_d_vec_2_1 = _mm512_and_si512(mask_vec, sr_p_vec_2_1);
            auto pre_d_vec_3_1 = _mm512_and_si512(mask_vec, sr_p_vec_3_1);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_2 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_3 = _mm512_loadu_si512((void*)(packed_loc + 48));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);
                    auto sl_p_vec_2 = _mm512_slli_epi32(p_vec_2, shift_amount);
                    auto sl_p_vec_3 = _mm512_slli_epi32(p_vec_3, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);
                    auto m_sl_p_vec_2 = _mm512_and_si512(mask_vec, sl_p_vec_2);
                    auto m_sl_p_vec_3 = _mm512_and_si512(mask_vec, sl_p_vec_3);

                    // add to result
                    pre_d_vec_0_1 = _mm512_or_si512(pre_d_vec_0_1, m_sl_p_vec_0);
                    pre_d_vec_1_1 = _mm512_or_si512(pre_d_vec_1_1, m_sl_p_vec_1);
                    pre_d_vec_2_1 = _mm512_or_si512(pre_d_vec_2_1, m_sl_p_vec_2);
                    pre_d_vec_3_1 = _mm512_or_si512(pre_d_vec_3_1, m_sl_p_vec_3);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto unsc_d_vec_0_1 = _mm512_slli_epi32(pre_d_vec_0_1, 16);
            auto unsc_d_vec_1_1 = _mm512_slli_epi32(pre_d_vec_1_1, 16);
            auto unsc_d_vec_2_1 = _mm512_slli_epi32(pre_d_vec_2_1, 16);
            auto unsc_d_vec_3_1 = _mm512_slli_epi32(pre_d_vec_3_1, 16);

            auto d_vec_0_1 = _mm512_mullo_epi16(unsc_d_vec_0_1, alpha_vec);
            auto d_vec_1_1 = _mm512_mullo_epi16(unsc_d_vec_1_1, alpha_vec);
            auto d_vec_2_1 = _mm512_mullo_epi16(unsc_d_vec_2_1, alpha_vec);
            auto d_vec_3_1 = _mm512_mullo_epi16(unsc_d_vec_3_1, alpha_vec);

            auto interleaved_0 = _mm512_or_si512(d_vec_0_0, d_vec_0_1);
            auto interleaved_1 = _mm512_or_si512(d_vec_1_0, d_vec_1_1);
            auto interleaved_2 = _mm512_or_si512(d_vec_2_0, d_vec_2_1);
            auto interleaved_3 = _mm512_or_si512(d_vec_3_0, d_vec_3_1);

            // store
            buf_type* store_loc = buf + (row >> 1)*buf_cols + (col << 1);
            _mm512_storeu_si512((void*)(store_loc), interleaved_0);
            _mm512_storeu_si512((void*)(store_loc + 32), interleaved_1);
            _mm512_storeu_si512((void*)(store_loc + 64), interleaved_2);
            _mm512_storeu_si512((void*)(store_loc + 96), interleaved_3);
        }

        if (row < row_num) {
            // init vec
            base_type* packed_loc = packed_data + col_loc.packed_num + col;
            p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
            p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));
            p_vec_2 = _mm512_loadu_si512((void*)(packed_loc + 32));
            p_vec_3 = _mm512_loadu_si512((void*)(packed_loc + 48));

            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);
            auto sr_p_vec_2 = _mm512_srli_epi32(p_vec_2, col_loc.bit_offset);
            auto sr_p_vec_3 = _mm512_srli_epi32(p_vec_3, col_loc.bit_offset);

            auto unsc_d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);
            auto unsc_d_vec_1 = _mm512_and_si512(mask_vec, sr_p_vec_1);
            auto unsc_d_vec_2 = _mm512_and_si512(mask_vec, sr_p_vec_2);
            auto unsc_d_vec_3 = _mm512_and_si512(mask_vec, sr_p_vec_3);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_2 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_3 = _mm512_loadu_si512((void*)(packed_loc + 48));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);
                    auto sl_p_vec_2 = _mm512_slli_epi32(p_vec_2, shift_amount);
                    auto sl_p_vec_3 = _mm512_slli_epi32(p_vec_3, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);
                    auto m_sl_p_vec_2 = _mm512_and_si512(mask_vec, sl_p_vec_2);
                    auto m_sl_p_vec_3 = _mm512_and_si512(mask_vec, sl_p_vec_3);

                    // add to result
                    unsc_d_vec_0 = _mm512_or_si512(unsc_d_vec_0, m_sl_p_vec_0);
                    unsc_d_vec_1 = _mm512_or_si512(unsc_d_vec_1, m_sl_p_vec_1);
                    unsc_d_vec_2 = _mm512_or_si512(unsc_d_vec_2, m_sl_p_vec_2);
                    unsc_d_vec_3 = _mm512_or_si512(unsc_d_vec_3, m_sl_p_vec_3);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto d_vec_0 = _mm512_mullo_epi16(unsc_d_vec_0, alpha_vec);
            auto d_vec_1 = _mm512_mullo_epi16(unsc_d_vec_1, alpha_vec);
            auto d_vec_2 = _mm512_mullo_epi16(unsc_d_vec_2, alpha_vec);
            auto d_vec_3 = _mm512_mullo_epi16(unsc_d_vec_3, alpha_vec);


            buf_type* store_loc = buf + (row >> 1)*buf_cols + (col << 1);
            _mm512_storeu_si512((void*)(store_loc), d_vec_0);
            _mm512_storeu_si512((void*)(store_loc + 32), d_vec_1);
            _mm512_storeu_si512((void*)(store_loc + 64), d_vec_2);
            _mm512_storeu_si512((void*)(store_loc + 96), d_vec_3);
        }
    }

    for (; col < col_num - 15; col+=16) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        base_type* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        int row;
        for (row = 0; row < row_num - 1; row+=2) {
            // vec 0
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto unsc_d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    unsc_d_vec_0_0 = _mm512_or_si512(unsc_d_vec_0_0, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto d_vec_0_0 = _mm512_mullo_epi16(unsc_d_vec_0_0, alpha_vec);

            // vec 1
            // shift and mask
            auto sr_p_vec_0_1 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto pre_d_vec_0_1 = _mm512_and_si512(mask_vec, sr_p_vec_0_1);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    pre_d_vec_0_1 = _mm512_or_si512(pre_d_vec_0_1, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto unsc_d_vec_0_1 = _mm512_slli_epi32(pre_d_vec_0_1, 16);
            auto d_vec_0_1 = _mm512_mullo_epi16(unsc_d_vec_0_1, alpha_vec);

            auto interleaved_0 = _mm512_or_si512(d_vec_0_0, d_vec_0_1);

            // store
            buf_type* store_loc = buf + (row >> 1)*buf_cols + (col << 1);
            _mm512_storeu_si512((void*)(store_loc), interleaved_0);
        }

        if (row < row_num) {
            // init vec
            base_type* packed_loc = packed_data + col_loc.packed_num + col;
            p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto unsc_d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    unsc_d_vec_0 = _mm512_or_si512(unsc_d_vec_0, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }


            auto d_vec_0 = _mm512_mullo_epi16(unsc_d_vec_0, alpha_vec);

            buf_type* store_loc = buf + (row >> 1)*buf_cols + (col << 1);
            _mm512_storeu_si512((void*)(store_loc), d_vec_0);
        }
    }

    if ((col < col_num) && (row_num & 1)) {
        memset((void*)(buf + (row_num >> 1)*buf_cols + (col << 1)), 0, (col_num - col) << 2);
    }

    // cleanup
    for (int row = 0; row < row_num; row++) {
        for (int col1 = col; col1 < col_num; col1++) {
            // get element
            buf[(row >> 1)*buf_cols + (col1 << 1) + (row & 1)] = mask & (buf_type)(packed_data[loc.packed_num + col1] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;

        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col1 = col; col1 < col_num; col1++) {
                    buf[(row >> 1)*buf_cols + (col1 << 1) + (row & 1)] |= mask & (packed_data[loc.packed_num + col1] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col1 = col; col1 < col_num; col1++) {
            buf[(row >> 1)*buf_cols + (col1 << 1) + (row & 1)] *= alpha;
        }
    }
}

// declare actual used specialization
template void scaled_madd_unpack<uint16_t, uint32_t>(uint16_t* buf, uint32_t* packed_data, uint16_t alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits);

/**
 * @brief Unpacks data in packed format into buffer in format for amx and scales each element by alpha
 * 
 * @tparam buf_type Type of destination buffer
 * @tparam base_type Basetype of APIB Mat object containing packed data
 * @param buf Pointer to destination buffer
 * @param packed_data Pointer to packed data
 * @param alpha 
 * @param loc Location of offset
 * @param row_num Number of rows to unpack
 * @param col_num Number of columns to unpack, assumes stays in same block
 * @param buf_cols Number of columns of destination buffer
 * @param num_bits Size of elements in bits
 */
template<typename buf_type, typename base_type>
void scaled_amx_unpack(buf_type* buf, base_type* packed_data, buf_type alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("matrix 32bit scale unpacking\n");
    uint32_t mask = (uint32_t)((1L << num_bits) - 1); 
    auto mask_vec = _mm512_set1_epi32(mask);
    auto alpha_vec = _mm512_set1_epi32(alpha);
    const int bits_type = 32;

    int col; 
    for (col = 0; col < col_num - 31; col+=32) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        base_type* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1 = _mm512_loadu_si512((void*)(packed_loc_init + 16));

        int row;
        for (row = 0; row < row_num - 3; row+=4) {
            // vec 0
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1_0 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);

            auto unsc_d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);
            auto unsc_d_vec_1_0 = _mm512_and_si512(mask_vec, sr_p_vec_1_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);

                    // add to result
                    unsc_d_vec_0_0 = _mm512_or_si512(unsc_d_vec_0_0, m_sl_p_vec_0);
                    unsc_d_vec_1_0 = _mm512_or_si512(unsc_d_vec_1_0, m_sl_p_vec_1);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto d_vec_0_0 = _mm512_mullo_epi32(unsc_d_vec_0_0, alpha_vec);
            auto d_vec_1_0 = _mm512_mullo_epi32(unsc_d_vec_1_0, alpha_vec);

            // vec 1
            // shift and mask
            auto sr_p_vec_0_1 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1_1 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);

            auto pre_d_vec_0_1 = _mm512_and_si512(mask_vec, sr_p_vec_0_1);
            auto pre_d_vec_1_1 = _mm512_and_si512(mask_vec, sr_p_vec_1_1);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);

                    // add to result
                    pre_d_vec_0_1 = _mm512_or_si512(pre_d_vec_0_1, m_sl_p_vec_0);
                    pre_d_vec_1_1 = _mm512_or_si512(pre_d_vec_1_1, m_sl_p_vec_1);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto sc_d_vec_0_1 = _mm512_mullo_epi32(pre_d_vec_0_1, alpha_vec);
            auto sc_d_vec_1_1 = _mm512_mullo_epi32(pre_d_vec_1_1, alpha_vec);

            auto d_vec_0_1 = _mm512_slli_epi32(sc_d_vec_0_1, 8);
            auto d_vec_1_1 = _mm512_slli_epi32(sc_d_vec_1_1, 8);

            // vec 2
            // shift and mask
            auto sr_p_vec_0_2 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1_2 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);

            auto pre_d_vec_0_2 = _mm512_and_si512(mask_vec, sr_p_vec_0_2);
            auto pre_d_vec_1_2 = _mm512_and_si512(mask_vec, sr_p_vec_1_2);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);

                    // add to result
                    pre_d_vec_0_2 = _mm512_or_si512(pre_d_vec_0_2, m_sl_p_vec_0);
                    pre_d_vec_1_2 = _mm512_or_si512(pre_d_vec_1_2, m_sl_p_vec_1);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto sc_d_vec_0_2 = _mm512_mullo_epi32(pre_d_vec_0_2, alpha_vec);
            auto sc_d_vec_1_2 = _mm512_mullo_epi32(pre_d_vec_1_2, alpha_vec);

            auto d_vec_0_2 = _mm512_slli_epi32(sc_d_vec_0_2, 16);
            auto d_vec_1_2 = _mm512_slli_epi32(sc_d_vec_1_2, 16);

            // vec 3
            // shift and mask
            auto sr_p_vec_0_3 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1_3 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);

            auto pre_d_vec_0_3 = _mm512_and_si512(mask_vec, sr_p_vec_0_3);
            auto pre_d_vec_1_3 = _mm512_and_si512(mask_vec, sr_p_vec_1_3);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);

                    // add to result
                    pre_d_vec_0_3 = _mm512_or_si512(pre_d_vec_0_3, m_sl_p_vec_0);
                    pre_d_vec_1_3 = _mm512_or_si512(pre_d_vec_1_3, m_sl_p_vec_1);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto sc_d_vec_0_3 = _mm512_mullo_epi32(pre_d_vec_0_3, alpha_vec);
            auto sc_d_vec_1_3 = _mm512_mullo_epi32(pre_d_vec_1_3, alpha_vec);

            auto d_vec_0_3 = _mm512_slli_epi32(sc_d_vec_0_3, 24);
            auto d_vec_1_3 = _mm512_slli_epi32(sc_d_vec_1_3, 24);

            auto comb_0_0 = _mm512_or_si512(d_vec_0_0, d_vec_0_1);
            auto comb_0_1 = _mm512_or_si512(d_vec_0_2, d_vec_0_3);

            auto comb_1_0 = _mm512_or_si512(d_vec_1_0, d_vec_1_1);
            auto comb_1_1 = _mm512_or_si512(d_vec_1_2, d_vec_1_3);

            auto interleaved_0 = _mm512_or_si512(comb_0_0, comb_0_1);
            auto interleaved_1 = _mm512_or_si512(comb_1_0, comb_1_1);

            // store
            buf_type* store_loc = buf + (row >> 2)*buf_cols + (col << 2);
            _mm512_storeu_si512((void*)(store_loc), interleaved_0);
            _mm512_storeu_si512((void*)(store_loc + 64), interleaved_1);
        }

        if (row < row_num) {
            // init vec
            base_type* packed_loc = packed_data + col_loc.packed_num + col;
            p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
            p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));
            auto res_vec_0 = _mm512_setzero_si512();
            auto res_vec_1 = _mm512_setzero_si512();
            int amnt = 0;

            for (;row < row_num; row++) {
                // shift and mask
                auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
                auto sr_p_vec_1 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);

                auto pre_d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);
                auto pre_d_vec_1 = _mm512_and_si512(mask_vec, sr_p_vec_1);

                col_loc.bit_offset += num_bits;

                // check for overlap
                if (col_loc.bit_offset >= bits_type) {
                    col_loc.packed_num += col_loc.eff_cols; 

                    // load next vec
                    base_type* packed_loc = packed_data + col_loc.packed_num + col;
                    p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                    p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));

                    if (col_loc.bit_offset > bits_type) {
                        // shift and mask
                        int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                        auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                        auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);

                        auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                        auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);

                        // add to result
                        pre_d_vec_0 = _mm512_or_si512(pre_d_vec_0, m_sl_p_vec_0);
                        pre_d_vec_1 = _mm512_or_si512(pre_d_vec_1, m_sl_p_vec_1);
                    }

                    col_loc.bit_offset &= (bits_type - 1);
                }

                auto sc_d_vec_0 = _mm512_mullo_epi32(pre_d_vec_0, alpha_vec);
                auto sc_d_vec_1 = _mm512_mullo_epi32(pre_d_vec_1, alpha_vec);

                auto d_vec_0 = _mm512_slli_epi32(sc_d_vec_0, amnt);
                auto d_vec_1 = _mm512_slli_epi32(sc_d_vec_1, amnt);

                res_vec_0 = _mm512_or_si512(d_vec_0, res_vec_0);
                res_vec_1 = _mm512_or_si512(d_vec_1, res_vec_1);

                amnt+=8;
            }

            buf_type* store_loc = buf + (row >> 2)*buf_cols + (col << 2);
            _mm512_storeu_si512((void*)(store_loc), res_vec_0);
            _mm512_storeu_si512((void*)(store_loc + 64), res_vec_1);
        }
    }

    for (; col < col_num - 15; col+=16) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        base_type* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        int row;
        for (row = 0; row < row_num - 3; row+=4) {
            // vec 0
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto unsc_d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    unsc_d_vec_0_0 = _mm512_or_si512(unsc_d_vec_0_0, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto d_vec_0_0 = _mm512_mullo_epi32(unsc_d_vec_0_0, alpha_vec);

            // vec 1
            // shift and mask
            auto sr_p_vec_0_1 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto pre_d_vec_0_1 = _mm512_and_si512(mask_vec, sr_p_vec_0_1);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    pre_d_vec_0_1 = _mm512_or_si512(pre_d_vec_0_1, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto sc_d_vec_0_1 = _mm512_mullo_epi32(pre_d_vec_0_1, alpha_vec);

            auto d_vec_0_1 = _mm512_slli_epi32(sc_d_vec_0_1, 8);

            // vec 2
            // shift and mask
            auto sr_p_vec_0_2 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto pre_d_vec_0_2 = _mm512_and_si512(mask_vec, sr_p_vec_0_2);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    pre_d_vec_0_2 = _mm512_or_si512(pre_d_vec_0_2, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto sc_d_vec_0_2 = _mm512_mullo_epi32(pre_d_vec_0_2, alpha_vec);

            auto d_vec_0_2 = _mm512_slli_epi32(sc_d_vec_0_2, 16);

            // vec 3
            // shift and mask
            auto sr_p_vec_0_3 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto pre_d_vec_0_3 = _mm512_and_si512(mask_vec, sr_p_vec_0_3);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                base_type* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    pre_d_vec_0_3 = _mm512_or_si512(pre_d_vec_0_3, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto sc_d_vec_0_3 = _mm512_mullo_epi32(pre_d_vec_0_3, alpha_vec);

            auto d_vec_0_3 = _mm512_slli_epi32(sc_d_vec_0_3, 24);

            auto comb_0_0 = _mm512_or_si512(d_vec_0_0, d_vec_0_1);
            auto comb_0_1 = _mm512_or_si512(d_vec_0_2, d_vec_0_3);

            auto interleaved_0 = _mm512_or_si512(comb_0_0, comb_0_1);

            // store
            buf_type* store_loc = buf + (row >> 2)*buf_cols + (col << 2);
            _mm512_storeu_si512((void*)(store_loc), interleaved_0);
        }

        if (row < row_num) {
            // init vec
            base_type* packed_loc = packed_data + col_loc.packed_num + col;
            p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
            auto res_vec_0 = _mm512_setzero_si512();
            int amnt = 0;

            for (;row < row_num; row++) {
                // shift and mask
                auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

                auto pre_d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);

                col_loc.bit_offset += num_bits;

                // check for overlap
                if (col_loc.bit_offset >= bits_type) {
                    col_loc.packed_num += col_loc.eff_cols; 

                    // load next vec
                    base_type* packed_loc = packed_data + col_loc.packed_num + col;
                    p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                    if (col_loc.bit_offset > bits_type) {
                        // shift and mask
                        int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                        auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                        auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                        // add to result
                        pre_d_vec_0 = _mm512_or_si512(pre_d_vec_0, m_sl_p_vec_0);
                    }

                    col_loc.bit_offset &= (bits_type - 1);
                }

                auto sc_d_vec_0 = _mm512_mullo_epi32(pre_d_vec_0, alpha_vec);

                auto d_vec_0 = _mm512_slli_epi32(sc_d_vec_0, amnt);

                res_vec_0 = _mm512_or_si512(d_vec_0, res_vec_0);

                amnt+=8;
            }

            buf_type* store_loc = buf + (row >> 2)*buf_cols + (col << 2);
            _mm512_storeu_si512((void*)(store_loc), res_vec_0);
        }
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        for (int col1 = col; col1 < col_num; col1++) {
            // get element
            buf[(row >> 2)*buf_cols + (col1 << 2) + (row & 3)] = mask & (buf_type)(packed_data[loc.packed_num + col1] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;

        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col1 = col; col1 < col_num; col1++) {
                    buf[(row >> 2)*buf_cols + (col1 << 2) + (row & 3)] |= mask & (packed_data[loc.packed_num + col1] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col1 = col; col1 < col_num; col1++) {
            buf[(row >> 2)*buf_cols + (col1 << 2) + (row & 3)] *= alpha;
        }
    }
}

template void scaled_amx_unpack<uint8_t, uint32_t>(uint8_t* buf, uint32_t* packed_data, uint8_t alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits);

/**
 * @brief Unpacks data in packed format into buffer in row-major format and scales each element by alpha
 * 
 * @tparam buf_type Type of destination buffer
 * @tparam base_type Basetype of APIB Mat object containing packed data
 * @param buf Pointer to destination buffer
 * @param packed_data Pointer to packed data
 * @param alpha 
 * @param loc Location of offset
 * @param row_num Number of rows to unpack
 * @param col_num Number of columns to unpack, assumes stays in same block
 * @param buf_cols Number of columns of destination buffer
 * @param num_bits Size of elements in bits
 */
template<typename buf_type, typename base_type>
void scaled_unpack(buf_type* buf, base_type* packed_data, buf_type alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("generic scale unpacking\n");
    buf_type mask = (buf_type)((1L << num_bits) - 1); 
    const int bits_type = 8 * sizeof(base_type);

    for (int row = 0; row < row_num; row++) {
        for (int col = 0; col < col_num; col++) {
            // get element
            buf[row*buf_cols + col] = mask & (buf_type)(packed_data[loc.packed_num + col] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;
        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col = 0; col < col_num; col++) {
                    buf[row*buf_cols + col] |= mask & (packed_data[loc.packed_num + col] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col = 0; col < col_num; col++) {
            buf[row*buf_cols + col] *= alpha;
        }
    }
}

// specialization to unpack and scale to 32 bits from basetype of 32 bits
template<>
void scaled_unpack(uint32_t* buf, uint32_t* packed_data, uint32_t alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("32bit scale unpacking\n");
    uint32_t mask = (uint32_t)((1L << num_bits) - 1); 
    auto mask_vec = _mm512_set1_epi32(mask);
    auto alpha_vec = _mm512_set1_epi32(alpha);
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
            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);
            auto sr_p_vec_2 = _mm512_srli_epi32(p_vec_2, col_loc.bit_offset);
            auto sr_p_vec_3 = _mm512_srli_epi32(p_vec_3, col_loc.bit_offset);
            auto sr_p_vec_4 = _mm512_srli_epi32(p_vec_4, col_loc.bit_offset);
            auto sr_p_vec_5 = _mm512_srli_epi32(p_vec_5, col_loc.bit_offset);
            auto sr_p_vec_6 = _mm512_srli_epi32(p_vec_6, col_loc.bit_offset);
            auto sr_p_vec_7 = _mm512_srli_epi32(p_vec_7, col_loc.bit_offset);

            auto d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);
            auto d_vec_1 = _mm512_and_si512(mask_vec, sr_p_vec_1);
            auto d_vec_2 = _mm512_and_si512(mask_vec, sr_p_vec_2);
            auto d_vec_3 = _mm512_and_si512(mask_vec, sr_p_vec_3);
            auto d_vec_4 = _mm512_and_si512(mask_vec, sr_p_vec_4);
            auto d_vec_5 = _mm512_and_si512(mask_vec, sr_p_vec_5);
            auto d_vec_6 = _mm512_and_si512(mask_vec, sr_p_vec_6);
            auto d_vec_7 = _mm512_and_si512(mask_vec, sr_p_vec_7);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_2 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_3 = _mm512_loadu_si512((void*)(packed_loc + 48));
                p_vec_4 = _mm512_loadu_si512((void*)(packed_loc + 64));
                p_vec_5 = _mm512_loadu_si512((void*)(packed_loc + 80));
                p_vec_6 = _mm512_loadu_si512((void*)(packed_loc + 96));
                p_vec_7 = _mm512_loadu_si512((void*)(packed_loc + 112));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);
                    auto sl_p_vec_2 = _mm512_slli_epi32(p_vec_2, shift_amount);
                    auto sl_p_vec_3 = _mm512_slli_epi32(p_vec_3, shift_amount);
                    auto sl_p_vec_4 = _mm512_slli_epi32(p_vec_4, shift_amount);
                    auto sl_p_vec_5 = _mm512_slli_epi32(p_vec_5, shift_amount);
                    auto sl_p_vec_6 = _mm512_slli_epi32(p_vec_6, shift_amount);
                    auto sl_p_vec_7 = _mm512_slli_epi32(p_vec_7, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);
                    auto m_sl_p_vec_2 = _mm512_and_si512(mask_vec, sl_p_vec_2);
                    auto m_sl_p_vec_3 = _mm512_and_si512(mask_vec, sl_p_vec_3);
                    auto m_sl_p_vec_4 = _mm512_and_si512(mask_vec, sl_p_vec_4);
                    auto m_sl_p_vec_5 = _mm512_and_si512(mask_vec, sl_p_vec_5);
                    auto m_sl_p_vec_6 = _mm512_and_si512(mask_vec, sl_p_vec_6);
                    auto m_sl_p_vec_7 = _mm512_and_si512(mask_vec, sl_p_vec_7);

                    // add to result
                    d_vec_0 = _mm512_or_si512(d_vec_0, m_sl_p_vec_0);
                    d_vec_1 = _mm512_or_si512(d_vec_1, m_sl_p_vec_1);
                    d_vec_2 = _mm512_or_si512(d_vec_2, m_sl_p_vec_2);
                    d_vec_3 = _mm512_or_si512(d_vec_3, m_sl_p_vec_3);
                    d_vec_4 = _mm512_or_si512(d_vec_4, m_sl_p_vec_4);
                    d_vec_5 = _mm512_or_si512(d_vec_5, m_sl_p_vec_5);
                    d_vec_6 = _mm512_or_si512(d_vec_6, m_sl_p_vec_6);
                    d_vec_7 = _mm512_or_si512(d_vec_7, m_sl_p_vec_7);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }
            // multiply
            auto scaled_vec_0 = _mm512_mullo_epi32(alpha_vec, d_vec_0);
            auto scaled_vec_1 = _mm512_mullo_epi32(alpha_vec, d_vec_1);
            auto scaled_vec_2 = _mm512_mullo_epi32(alpha_vec, d_vec_2);
            auto scaled_vec_3 = _mm512_mullo_epi32(alpha_vec, d_vec_3);
            auto scaled_vec_4 = _mm512_mullo_epi32(alpha_vec, d_vec_4);
            auto scaled_vec_5 = _mm512_mullo_epi32(alpha_vec, d_vec_5);
            auto scaled_vec_6 = _mm512_mullo_epi32(alpha_vec, d_vec_6);
            auto scaled_vec_7 = _mm512_mullo_epi32(alpha_vec, d_vec_7);

            // store
            uint32_t* store_loc = buf + row*buf_cols + col;
            _mm512_storeu_si512((void*)(store_loc), scaled_vec_0);
            _mm512_storeu_si512((void*)(store_loc + 16), scaled_vec_1);
            _mm512_storeu_si512((void*)(store_loc + 32), scaled_vec_2);
            _mm512_storeu_si512((void*)(store_loc + 48), scaled_vec_3);
            _mm512_storeu_si512((void*)(store_loc + 64), scaled_vec_4);
            _mm512_storeu_si512((void*)(store_loc + 80), scaled_vec_5);
            _mm512_storeu_si512((void*)(store_loc + 96), scaled_vec_6);
            _mm512_storeu_si512((void*)(store_loc + 112), scaled_vec_7);
        }
    }

    // first cleanup
    for (; col < col_num - 15; col+=16) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        for (int row = 0; row < row_num; row++) {
            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    d_vec_0 = _mm512_or_si512(d_vec_0, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }
            // multiply
            auto scaled_vec_0 = _mm512_mullo_epi32(alpha_vec, d_vec_0);

            // store
            uint32_t* store_loc = buf + row*buf_cols + col;
            _mm512_storeu_si512((void*)(store_loc), scaled_vec_0);
        }
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        for (int col1 = col; col1 < col_num; col1++) {
            // get element
            buf[row*buf_cols + col1] = mask & (uint32_t)(packed_data[loc.packed_num + col1] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;

        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col1 = col; col1 < col_num; col1++) {
                    buf[row*buf_cols + col1] |= mask & (packed_data[loc.packed_num + col1] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col1 = col; col1 < col_num; col1++) {
            buf[row*buf_cols + col1] *= alpha;
        }
    }
}

// specialization to unpack and scale to 16 bits from basetype of 32 bits
template<>
void scaled_unpack(uint16_t* buf, uint32_t* packed_data, uint16_t alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("16bit scale unpacking\n");
    uint32_t mask = (uint32_t)((1L << num_bits) - 1); 
    auto mask_vec = _mm512_set1_epi32(mask);
    auto alpha_vec = _mm512_set1_epi16(alpha);
    const int bits_type = 32;

    auto permute_vec = _mm512_set_epi64(7, 5, 3, 1, 6, 4, 2, 0);

    int col; 
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
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0_0, col_loc.bit_offset);
            auto sr_p_vec_1_0 = _mm512_srli_epi32(p_vec_1_0, col_loc.bit_offset);

            auto sr_p_vec_0_1 = _mm512_srli_epi32(p_vec_0_1, col_loc.bit_offset);
            auto sr_p_vec_1_1 = _mm512_srli_epi32(p_vec_1_1, col_loc.bit_offset);

            auto sr_p_vec_0_2 = _mm512_srli_epi32(p_vec_0_2, col_loc.bit_offset);
            auto sr_p_vec_1_2 = _mm512_srli_epi32(p_vec_1_2, col_loc.bit_offset);

            auto sr_p_vec_0_3 = _mm512_srli_epi32(p_vec_0_3, col_loc.bit_offset);
            auto sr_p_vec_1_3 = _mm512_srli_epi32(p_vec_1_3, col_loc.bit_offset);

            auto d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);
            auto d_vec_1_0 = _mm512_and_si512(mask_vec, sr_p_vec_1_0);

            auto d_vec_0_1 = _mm512_and_si512(mask_vec, sr_p_vec_0_1);
            auto d_vec_1_1 = _mm512_and_si512(mask_vec, sr_p_vec_1_1);

            auto d_vec_0_2 = _mm512_and_si512(mask_vec, sr_p_vec_0_2);
            auto d_vec_1_2 = _mm512_and_si512(mask_vec, sr_p_vec_1_2);

            auto d_vec_0_3 = _mm512_and_si512(mask_vec, sr_p_vec_0_3);
            auto d_vec_1_3 = _mm512_and_si512(mask_vec, sr_p_vec_1_3);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc + 16));

                p_vec_0_1 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_1_1 = _mm512_loadu_si512((void*)(packed_loc + 48));

                p_vec_0_2 = _mm512_loadu_si512((void*)(packed_loc + 64));
                p_vec_1_2 = _mm512_loadu_si512((void*)(packed_loc + 80));

                p_vec_0_3 = _mm512_loadu_si512((void*)(packed_loc + 96));
                p_vec_1_3 = _mm512_loadu_si512((void*)(packed_loc + 112));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0_0 = _mm512_slli_epi32(p_vec_0_0, shift_amount);
                    auto sl_p_vec_1_0 = _mm512_slli_epi32(p_vec_1_0, shift_amount);

                    auto sl_p_vec_0_1 = _mm512_slli_epi32(p_vec_0_1, shift_amount);
                    auto sl_p_vec_1_1 = _mm512_slli_epi32(p_vec_1_1, shift_amount);

                    auto sl_p_vec_0_2 = _mm512_slli_epi32(p_vec_0_2, shift_amount);
                    auto sl_p_vec_1_2 = _mm512_slli_epi32(p_vec_1_2, shift_amount);

                    auto sl_p_vec_0_3 = _mm512_slli_epi32(p_vec_0_3, shift_amount);
                    auto sl_p_vec_1_3 = _mm512_slli_epi32(p_vec_1_3, shift_amount);

                    auto m_sl_p_vec_0_0 = _mm512_and_si512(mask_vec, sl_p_vec_0_0);
                    auto m_sl_p_vec_1_0 = _mm512_and_si512(mask_vec, sl_p_vec_1_0);

                    auto m_sl_p_vec_0_1 = _mm512_and_si512(mask_vec, sl_p_vec_0_1);
                    auto m_sl_p_vec_1_1 = _mm512_and_si512(mask_vec, sl_p_vec_1_1);

                    auto m_sl_p_vec_0_2 = _mm512_and_si512(mask_vec, sl_p_vec_0_2);
                    auto m_sl_p_vec_1_2 = _mm512_and_si512(mask_vec, sl_p_vec_1_2);

                    auto m_sl_p_vec_0_3 = _mm512_and_si512(mask_vec, sl_p_vec_0_3);
                    auto m_sl_p_vec_1_3 = _mm512_and_si512(mask_vec, sl_p_vec_1_3);

                    // add to result
                    d_vec_0_0 = _mm512_or_si512(d_vec_0_0, m_sl_p_vec_0_0);
                    d_vec_1_0 = _mm512_or_si512(d_vec_1_0, m_sl_p_vec_1_0);

                    d_vec_0_1 = _mm512_or_si512(d_vec_0_1, m_sl_p_vec_0_1);
                    d_vec_1_1 = _mm512_or_si512(d_vec_1_1, m_sl_p_vec_1_1);

                    d_vec_0_2 = _mm512_or_si512(d_vec_0_2, m_sl_p_vec_0_2);
                    d_vec_1_2 = _mm512_or_si512(d_vec_1_2, m_sl_p_vec_1_2);

                    d_vec_0_3 = _mm512_or_si512(d_vec_0_3, m_sl_p_vec_0_3);
                    d_vec_1_3 = _mm512_or_si512(d_vec_1_3, m_sl_p_vec_1_3);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            // store
            uint16_t* store_loc = buf + row*buf_cols + col;
            auto packed_0 = _mm512_packus_epi32(d_vec_0_0, d_vec_1_0);
            auto packed_1 = _mm512_packus_epi32(d_vec_0_1, d_vec_1_1);
            auto packed_2 = _mm512_packus_epi32(d_vec_0_2, d_vec_1_2);
            auto packed_3 = _mm512_packus_epi32(d_vec_0_3, d_vec_1_3);

            auto aligned_0 = _mm512_permutexvar_epi64(permute_vec, packed_0);
            auto aligned_1 = _mm512_permutexvar_epi64(permute_vec, packed_1);
            auto aligned_2 = _mm512_permutexvar_epi64(permute_vec, packed_2);
            auto aligned_3 = _mm512_permutexvar_epi64(permute_vec, packed_3);

            auto scaled_0 = _mm512_mullo_epi16(alpha_vec, aligned_0);
            auto scaled_1 = _mm512_mullo_epi16(alpha_vec, aligned_1);
            auto scaled_2 = _mm512_mullo_epi16(alpha_vec, aligned_2);
            auto scaled_3 = _mm512_mullo_epi16(alpha_vec, aligned_3);

            _mm512_storeu_si512((void*)(store_loc), scaled_0);
            _mm512_storeu_si512((void*)(store_loc + 32), scaled_1);
            _mm512_storeu_si512((void*)(store_loc + 64), scaled_2);
            _mm512_storeu_si512((void*)(store_loc + 96), scaled_3);

        }
    }

    // first cleanup
    for (; col < col_num - 31; col+=32) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc_init));
        auto p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc_init + 16));

        for (int row = 0; row < row_num; row++) {
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0_0, col_loc.bit_offset);
            auto sr_p_vec_1_0 = _mm512_srli_epi32(p_vec_1_0, col_loc.bit_offset);

            auto d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);
            auto d_vec_1_0 = _mm512_and_si512(mask_vec, sr_p_vec_1_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc + 16));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0_0 = _mm512_slli_epi32(p_vec_0_0, shift_amount);
                    auto sl_p_vec_1_0 = _mm512_slli_epi32(p_vec_1_0, shift_amount);

                    auto m_sl_p_vec_0_0 = _mm512_and_si512(mask_vec, sl_p_vec_0_0);
                    auto m_sl_p_vec_1_0 = _mm512_and_si512(mask_vec, sl_p_vec_1_0);

                    // add to result
                    d_vec_0_0 = _mm512_or_si512(d_vec_0_0, m_sl_p_vec_0_0);
                    d_vec_1_0 = _mm512_or_si512(d_vec_1_0, m_sl_p_vec_1_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            // store
            uint16_t* store_loc = buf + row*buf_cols + col;
            auto packed_0 = _mm512_packus_epi32(d_vec_0_0, d_vec_1_0);
            auto aligned_0 = _mm512_permutexvar_epi64(permute_vec, packed_0);
            auto scaled_0 = _mm512_mullo_epi16(alpha_vec, aligned_0);

            _mm512_storeu_si512((void*)(store_loc), scaled_0);

        }
    }

    // cleanup
    for (int row = 0; row < row_num; row++) {
        for (int col1 = col; col1 < col_num; col1++) {
            // get element
            buf[row*buf_cols + col1] = mask & (uint16_t)(packed_data[loc.packed_num + col1] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;

        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col1 = col; col1 < col_num; col1++) {
                    buf[row*buf_cols + col1] |= mask & (packed_data[loc.packed_num + col1] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col1 = col; col1 < col_num; col1++) {
            buf[row*buf_cols + col1] *= alpha;
        }
    }
}

// specialization to unpack and scale to 8 bits from basetype of 32 bits
template<>
void scaled_unpack<uint8_t, uint32_t>(uint8_t* buf, uint32_t* packed_data, uint8_t alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("8bit scale unpacking\n");
    uint32_t mask = (uint32_t)((1L << num_bits) - 1); 
    auto mask_vec = _mm512_set1_epi32(mask);
    auto permute_vec = _mm512_set_epi32(15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0);
    auto alpha_vec = _mm512_set1_epi16(alpha);
    const int bits_type = 32;

    int col; 
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
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0_0, col_loc.bit_offset);
            auto sr_p_vec_1_0 = _mm512_srli_epi32(p_vec_1_0, col_loc.bit_offset);
            auto sr_p_vec_2_0 = _mm512_srli_epi32(p_vec_2_0, col_loc.bit_offset);
            auto sr_p_vec_3_0 = _mm512_srli_epi32(p_vec_3_0, col_loc.bit_offset);

            auto sr_p_vec_0_1 = _mm512_srli_epi32(p_vec_0_1, col_loc.bit_offset);
            auto sr_p_vec_1_1 = _mm512_srli_epi32(p_vec_1_1, col_loc.bit_offset);
            auto sr_p_vec_2_1 = _mm512_srli_epi32(p_vec_2_1, col_loc.bit_offset);
            auto sr_p_vec_3_1 = _mm512_srli_epi32(p_vec_3_1, col_loc.bit_offset);

            auto d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);
            auto d_vec_1_0 = _mm512_and_si512(mask_vec, sr_p_vec_1_0);
            auto d_vec_2_0 = _mm512_and_si512(mask_vec, sr_p_vec_2_0);
            auto d_vec_3_0 = _mm512_and_si512(mask_vec, sr_p_vec_3_0);

            auto d_vec_0_1 = _mm512_and_si512(mask_vec, sr_p_vec_0_1);
            auto d_vec_1_1 = _mm512_and_si512(mask_vec, sr_p_vec_1_1);
            auto d_vec_2_1 = _mm512_and_si512(mask_vec, sr_p_vec_2_1);
            auto d_vec_3_1 = _mm512_and_si512(mask_vec, sr_p_vec_3_1);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_2_0 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_3_0 = _mm512_loadu_si512((void*)(packed_loc + 48));

                p_vec_0_1 = _mm512_loadu_si512((void*)(packed_loc + 64));
                p_vec_1_1 = _mm512_loadu_si512((void*)(packed_loc + 80));
                p_vec_2_1 = _mm512_loadu_si512((void*)(packed_loc + 96));
                p_vec_3_1 = _mm512_loadu_si512((void*)(packed_loc + 112));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0_0 = _mm512_slli_epi32(p_vec_0_0, shift_amount);
                    auto sl_p_vec_1_0 = _mm512_slli_epi32(p_vec_1_0, shift_amount);
                    auto sl_p_vec_2_0 = _mm512_slli_epi32(p_vec_2_0, shift_amount);
                    auto sl_p_vec_3_0 = _mm512_slli_epi32(p_vec_3_0, shift_amount);

                    auto sl_p_vec_0_1 = _mm512_slli_epi32(p_vec_0_1, shift_amount);
                    auto sl_p_vec_1_1 = _mm512_slli_epi32(p_vec_1_1, shift_amount);
                    auto sl_p_vec_2_1 = _mm512_slli_epi32(p_vec_2_1, shift_amount);
                    auto sl_p_vec_3_1 = _mm512_slli_epi32(p_vec_3_1, shift_amount);

                    auto m_sl_p_vec_0_0 = _mm512_and_si512(mask_vec, sl_p_vec_0_0);
                    auto m_sl_p_vec_1_0 = _mm512_and_si512(mask_vec, sl_p_vec_1_0);
                    auto m_sl_p_vec_2_0 = _mm512_and_si512(mask_vec, sl_p_vec_2_0);
                    auto m_sl_p_vec_3_0 = _mm512_and_si512(mask_vec, sl_p_vec_3_0);

                    auto m_sl_p_vec_0_1 = _mm512_and_si512(mask_vec, sl_p_vec_0_1);
                    auto m_sl_p_vec_1_1 = _mm512_and_si512(mask_vec, sl_p_vec_1_1);
                    auto m_sl_p_vec_2_1 = _mm512_and_si512(mask_vec, sl_p_vec_2_1);
                    auto m_sl_p_vec_3_1 = _mm512_and_si512(mask_vec, sl_p_vec_3_1);

                    // add to result
                    d_vec_0_0 = _mm512_or_si512(d_vec_0_0, m_sl_p_vec_0_0);
                    d_vec_1_0 = _mm512_or_si512(d_vec_1_0, m_sl_p_vec_1_0);
                    d_vec_2_0 = _mm512_or_si512(d_vec_2_0, m_sl_p_vec_2_0);
                    d_vec_3_0 = _mm512_or_si512(d_vec_3_0, m_sl_p_vec_3_0);

                    d_vec_0_1 = _mm512_or_si512(d_vec_0_1, m_sl_p_vec_0_1);
                    d_vec_1_1 = _mm512_or_si512(d_vec_1_1, m_sl_p_vec_1_1);
                    d_vec_2_1 = _mm512_or_si512(d_vec_2_1, m_sl_p_vec_2_1);
                    d_vec_3_1 = _mm512_or_si512(d_vec_3_1, m_sl_p_vec_3_1);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            // store
            uint8_t* store_loc = buf + row*buf_cols + col;
            auto packed_16_0_0 = _mm512_packus_epi32(d_vec_0_0, d_vec_1_0);
            auto packed_16_1_0 = _mm512_packus_epi32(d_vec_2_0, d_vec_3_0);

            auto packed_16_0_1 = _mm512_packus_epi32(d_vec_0_1, d_vec_1_1);
            auto packed_16_1_1 = _mm512_packus_epi32(d_vec_2_1, d_vec_3_1);

            auto scaled_16_0_0 = _mm512_mullo_epi16(alpha_vec, packed_16_0_0);
            auto scaled_16_1_0 = _mm512_mullo_epi16(alpha_vec, packed_16_1_0);

            auto scaled_16_0_1 = _mm512_mullo_epi16(alpha_vec, packed_16_0_1);
            auto scaled_16_1_1 = _mm512_mullo_epi16(alpha_vec, packed_16_1_1);

            auto packed_0 = _mm512_packus_epi16(scaled_16_0_0, scaled_16_1_0);
            auto packed_1 = _mm512_packus_epi16(scaled_16_0_1, scaled_16_1_1);

            auto aligned_0 = _mm512_permutexvar_epi32(permute_vec, packed_0);
            auto aligned_1 = _mm512_permutexvar_epi32(permute_vec, packed_1);

            _mm512_storeu_si512((void*)(store_loc), aligned_0);
            _mm512_storeu_si512((void*)(store_loc + 64), aligned_1);

        }
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
            // shift and mask
            auto sr_p_vec_0_0 = _mm512_srli_epi32(p_vec_0_0, col_loc.bit_offset);
            auto sr_p_vec_1_0 = _mm512_srli_epi32(p_vec_1_0, col_loc.bit_offset);
            auto sr_p_vec_2_0 = _mm512_srli_epi32(p_vec_2_0, col_loc.bit_offset);
            auto sr_p_vec_3_0 = _mm512_srli_epi32(p_vec_3_0, col_loc.bit_offset);

            auto d_vec_0_0 = _mm512_and_si512(mask_vec, sr_p_vec_0_0);
            auto d_vec_1_0 = _mm512_and_si512(mask_vec, sr_p_vec_1_0);
            auto d_vec_2_0 = _mm512_and_si512(mask_vec, sr_p_vec_2_0);
            auto d_vec_3_0 = _mm512_and_si512(mask_vec, sr_p_vec_3_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1_0 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_2_0 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_3_0 = _mm512_loadu_si512((void*)(packed_loc + 48));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0_0 = _mm512_slli_epi32(p_vec_0_0, shift_amount);
                    auto sl_p_vec_1_0 = _mm512_slli_epi32(p_vec_1_0, shift_amount);
                    auto sl_p_vec_2_0 = _mm512_slli_epi32(p_vec_2_0, shift_amount);
                    auto sl_p_vec_3_0 = _mm512_slli_epi32(p_vec_3_0, shift_amount);

                    auto m_sl_p_vec_0_0 = _mm512_and_si512(mask_vec, sl_p_vec_0_0);
                    auto m_sl_p_vec_1_0 = _mm512_and_si512(mask_vec, sl_p_vec_1_0);
                    auto m_sl_p_vec_2_0 = _mm512_and_si512(mask_vec, sl_p_vec_2_0);
                    auto m_sl_p_vec_3_0 = _mm512_and_si512(mask_vec, sl_p_vec_3_0);

                    // add to result
                    d_vec_0_0 = _mm512_or_si512(d_vec_0_0, m_sl_p_vec_0_0);
                    d_vec_1_0 = _mm512_or_si512(d_vec_1_0, m_sl_p_vec_1_0);
                    d_vec_2_0 = _mm512_or_si512(d_vec_2_0, m_sl_p_vec_2_0);
                    d_vec_3_0 = _mm512_or_si512(d_vec_3_0, m_sl_p_vec_3_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            // store
            uint8_t* store_loc = buf + row*buf_cols + col;
            auto packed_16_0_0 = _mm512_packus_epi32(d_vec_0_0, d_vec_1_0);
            auto packed_16_1_0 = _mm512_packus_epi32(d_vec_2_0, d_vec_3_0);

            auto scaled_16_0_0 = _mm512_mullo_epi16(alpha_vec, packed_16_0_0);
            auto scaled_16_1_0 = _mm512_mullo_epi16(alpha_vec, packed_16_1_0);

            auto packed_0 = _mm512_packus_epi16(scaled_16_0_0, scaled_16_1_0);

            auto aligned_0 = _mm512_permutexvar_epi32(permute_vec, packed_0);

            _mm512_storeu_si512((void*)(store_loc), aligned_0);

        }
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        for (int col1 = col; col1 < col_num; col1++) {
            // get element
            buf[row*buf_cols + col1] = mask & (uint32_t)(packed_data[loc.packed_num + col1] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;

        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col1 = col; col1 < col_num; col1++) {
                    buf[row*buf_cols + col1] |= mask & (packed_data[loc.packed_num + col1] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col1 = col; col1 < col_num; col1++) {
            buf[row*buf_cols + col1] *= alpha;
        }
    }
}

// specialization to unpack and scale to 64 bits from basetype of 32 bits
template<>
void scaled_unpack<uint64_t, uint32_t>(uint64_t* buf, uint32_t* packed_data, uint64_t alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("64bit scale unpacking\n");
    uint32_t mask = (uint32_t)((1L << num_bits) - 1); 
    auto mask_vec = _mm512_set1_epi32(mask);
    auto permute_vec = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    auto alpha_vec = _mm512_set1_epi64(alpha);
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
            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1 = _mm512_srli_epi32(p_vec_1, col_loc.bit_offset);
            auto sr_p_vec_2 = _mm512_srli_epi32(p_vec_2, col_loc.bit_offset);
            auto sr_p_vec_3 = _mm512_srli_epi32(p_vec_3, col_loc.bit_offset);

            auto d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);
            auto d_vec_1 = _mm512_and_si512(mask_vec, sr_p_vec_1);
            auto d_vec_2 = _mm512_and_si512(mask_vec, sr_p_vec_2);
            auto d_vec_3 = _mm512_and_si512(mask_vec, sr_p_vec_3);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_2 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_3 = _mm512_loadu_si512((void*)(packed_loc + 48));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi32(p_vec_1, shift_amount);
                    auto sl_p_vec_2 = _mm512_slli_epi32(p_vec_2, shift_amount);
                    auto sl_p_vec_3 = _mm512_slli_epi32(p_vec_3, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);
                    auto m_sl_p_vec_2 = _mm512_and_si512(mask_vec, sl_p_vec_2);
                    auto m_sl_p_vec_3 = _mm512_and_si512(mask_vec, sl_p_vec_3);

                    // add to result
                    d_vec_0 = _mm512_or_si512(d_vec_0, m_sl_p_vec_0);
                    d_vec_1 = _mm512_or_si512(d_vec_1, m_sl_p_vec_1);
                    d_vec_2 = _mm512_or_si512(d_vec_2, m_sl_p_vec_2);
                    d_vec_3 = _mm512_or_si512(d_vec_3, m_sl_p_vec_3);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto d_vec_0_0 = _mm512_castsi512_si256(d_vec_0);
            auto d_vec_1_0 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_0));

            auto d_vec_0_1 = _mm512_castsi512_si256(d_vec_1);
            auto d_vec_1_1 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_1));

            auto d_vec_0_2 = _mm512_castsi512_si256(d_vec_2);
            auto d_vec_1_2 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_2));

            auto d_vec_0_3 = _mm512_castsi512_si256(d_vec_3);
            auto d_vec_1_3 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_3));

            auto d_vec_0_0_64 = _mm512_cvtepu32_epi64(d_vec_0_0);
            auto d_vec_1_0_64 = _mm512_cvtepu32_epi64(d_vec_1_0);

            auto d_vec_0_1_64 = _mm512_cvtepu32_epi64(d_vec_0_1);
            auto d_vec_1_1_64 = _mm512_cvtepu32_epi64(d_vec_1_1);

            auto d_vec_0_2_64 = _mm512_cvtepu32_epi64(d_vec_0_2);
            auto d_vec_1_2_64 = _mm512_cvtepu32_epi64(d_vec_1_2);

            auto d_vec_0_3_64 = _mm512_cvtepu32_epi64(d_vec_0_3);
            auto d_vec_1_3_64 = _mm512_cvtepu32_epi64(d_vec_1_3);

            auto d_vec_0_0_scaled = _mm512_mullo_epi64(d_vec_0_0_64, alpha_vec);
            auto d_vec_1_0_scaled = _mm512_mullo_epi64(d_vec_1_0_64, alpha_vec);

            auto d_vec_0_1_scaled = _mm512_mullo_epi64(d_vec_0_1_64, alpha_vec);
            auto d_vec_1_1_scaled = _mm512_mullo_epi64(d_vec_1_1_64, alpha_vec);

            auto d_vec_0_2_scaled = _mm512_mullo_epi64(d_vec_0_2_64, alpha_vec);
            auto d_vec_1_2_scaled = _mm512_mullo_epi64(d_vec_1_2_64, alpha_vec);

            auto d_vec_0_3_scaled = _mm512_mullo_epi64(d_vec_0_3_64, alpha_vec);
            auto d_vec_1_3_scaled = _mm512_mullo_epi64(d_vec_1_3_64, alpha_vec);

            // store
            uint64_t* store_loc = buf + row*buf_cols + col;
            _mm512_storeu_si512((void*)(store_loc), d_vec_0_0_scaled);
            _mm512_storeu_si512((void*)(store_loc + 8), d_vec_1_0_scaled);
            _mm512_storeu_si512((void*)(store_loc + 16), d_vec_0_1_scaled);
            _mm512_storeu_si512((void*)(store_loc + 24), d_vec_1_1_scaled);
            _mm512_storeu_si512((void*)(store_loc + 32), d_vec_0_2_scaled);
            _mm512_storeu_si512((void*)(store_loc + 40), d_vec_1_2_scaled);
            _mm512_storeu_si512((void*)(store_loc + 48), d_vec_0_3_scaled);
            _mm512_storeu_si512((void*)(store_loc + 56), d_vec_1_3_scaled);
        }
    }

    // first cleanup
    for (; col < col_num - 15; col+=16) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        uint32_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        for (int row = 0; row < row_num; row++) {
            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi32(p_vec_0, col_loc.bit_offset);

            auto d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint32_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi32(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    d_vec_0 = _mm512_or_si512(d_vec_0, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto d_vec_0_0 = _mm512_castsi512_si256(d_vec_0);
            auto d_vec_1_0 = _mm512_castsi512_si256(_mm512_permutexvar_epi64(permute_vec, d_vec_0));

            auto d_vec_0_0_64 = _mm512_cvtepu32_epi64(d_vec_0_0);
            auto d_vec_1_0_64 = _mm512_cvtepu32_epi64(d_vec_1_0);

            auto d_vec_0_0_scaled = _mm512_mullo_epi64(d_vec_0_0_64, alpha_vec);
            auto d_vec_1_0_scaled = _mm512_mullo_epi64(d_vec_1_0_64, alpha_vec);

            // store
            uint64_t* store_loc = buf + row*buf_cols + col;

            _mm512_storeu_si512((void*)(store_loc), d_vec_0_0_scaled);
            _mm512_storeu_si512((void*)(store_loc + 8), d_vec_1_0_scaled);
        }
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        for (int col1 = col; col1 < col_num; col1++) {
            // get element
            buf[row*buf_cols + col1] = mask & (uint32_t)(packed_data[loc.packed_num + col1] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;

        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col1 = col; col1 < col_num; col1++) {
                    buf[row*buf_cols + col1] |= mask & (uint32_t)(packed_data[loc.packed_num + col1] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col1 = col; col1 < col_num; col1++) {
            buf[row*buf_cols + col1] *= alpha;
        }
    }
}

// specialization to unpack and scale to 64 bits from basetype of 64 bits
template<>
void scaled_unpack<uint64_t, uint64_t>(uint64_t* buf, uint64_t* packed_data, uint64_t alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits) {
    if (CHECK_TEMPLATE) printf("64bit to 64bit scale unpacking\n");
    uint64_t mask;
    if (num_bits >= 64) {
        mask = (uint64_t)(-1);
    } else if (num_bits >= 63) {
        mask = (uint64_t)(~(1L << num_bits));
    } else {
        mask = (uint64_t)((1L << num_bits) - 1); 
    }

    auto mask_vec = _mm512_set1_epi64(mask);
    const int bits_type = 64;
    auto alpha_vec = _mm512_set1_epi64(alpha);

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
            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi64(p_vec_0, col_loc.bit_offset);
            auto sr_p_vec_1 = _mm512_srli_epi64(p_vec_1, col_loc.bit_offset);
            auto sr_p_vec_2 = _mm512_srli_epi64(p_vec_2, col_loc.bit_offset);
            auto sr_p_vec_3 = _mm512_srli_epi64(p_vec_3, col_loc.bit_offset);
            auto sr_p_vec_4 = _mm512_srli_epi64(p_vec_4, col_loc.bit_offset);
            auto sr_p_vec_5 = _mm512_srli_epi64(p_vec_5, col_loc.bit_offset);
            auto sr_p_vec_6 = _mm512_srli_epi64(p_vec_6, col_loc.bit_offset);
            auto sr_p_vec_7 = _mm512_srli_epi64(p_vec_7, col_loc.bit_offset);

            auto d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);
            auto d_vec_1 = _mm512_and_si512(mask_vec, sr_p_vec_1);
            auto d_vec_2 = _mm512_and_si512(mask_vec, sr_p_vec_2);
            auto d_vec_3 = _mm512_and_si512(mask_vec, sr_p_vec_3);
            auto d_vec_4 = _mm512_and_si512(mask_vec, sr_p_vec_4);
            auto d_vec_5 = _mm512_and_si512(mask_vec, sr_p_vec_5);
            auto d_vec_6 = _mm512_and_si512(mask_vec, sr_p_vec_6);
            auto d_vec_7 = _mm512_and_si512(mask_vec, sr_p_vec_7);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint64_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));
                p_vec_1 = _mm512_loadu_si512((void*)(packed_loc + 8));
                p_vec_2 = _mm512_loadu_si512((void*)(packed_loc + 16));
                p_vec_3 = _mm512_loadu_si512((void*)(packed_loc + 24));
                p_vec_4 = _mm512_loadu_si512((void*)(packed_loc + 32));
                p_vec_5 = _mm512_loadu_si512((void*)(packed_loc + 40));
                p_vec_6 = _mm512_loadu_si512((void*)(packed_loc + 48));
                p_vec_7 = _mm512_loadu_si512((void*)(packed_loc + 56));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi64(p_vec_0, shift_amount);
                    auto sl_p_vec_1 = _mm512_slli_epi64(p_vec_1, shift_amount);
                    auto sl_p_vec_2 = _mm512_slli_epi64(p_vec_2, shift_amount);
                    auto sl_p_vec_3 = _mm512_slli_epi64(p_vec_3, shift_amount);
                    auto sl_p_vec_4 = _mm512_slli_epi64(p_vec_4, shift_amount);
                    auto sl_p_vec_5 = _mm512_slli_epi64(p_vec_5, shift_amount);
                    auto sl_p_vec_6 = _mm512_slli_epi64(p_vec_6, shift_amount);
                    auto sl_p_vec_7 = _mm512_slli_epi64(p_vec_7, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);
                    auto m_sl_p_vec_1 = _mm512_and_si512(mask_vec, sl_p_vec_1);
                    auto m_sl_p_vec_2 = _mm512_and_si512(mask_vec, sl_p_vec_2);
                    auto m_sl_p_vec_3 = _mm512_and_si512(mask_vec, sl_p_vec_3);
                    auto m_sl_p_vec_4 = _mm512_and_si512(mask_vec, sl_p_vec_4);
                    auto m_sl_p_vec_5 = _mm512_and_si512(mask_vec, sl_p_vec_5);
                    auto m_sl_p_vec_6 = _mm512_and_si512(mask_vec, sl_p_vec_6);
                    auto m_sl_p_vec_7 = _mm512_and_si512(mask_vec, sl_p_vec_7);

                    // add to result
                    d_vec_0 = _mm512_or_si512(d_vec_0, m_sl_p_vec_0);
                    d_vec_1 = _mm512_or_si512(d_vec_1, m_sl_p_vec_1);
                    d_vec_2 = _mm512_or_si512(d_vec_2, m_sl_p_vec_2);
                    d_vec_3 = _mm512_or_si512(d_vec_3, m_sl_p_vec_3);
                    d_vec_4 = _mm512_or_si512(d_vec_4, m_sl_p_vec_4);
                    d_vec_5 = _mm512_or_si512(d_vec_5, m_sl_p_vec_5);
                    d_vec_6 = _mm512_or_si512(d_vec_6, m_sl_p_vec_6);
                    d_vec_7 = _mm512_or_si512(d_vec_7, m_sl_p_vec_7);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto scaled_d_vec_0 = _mm512_mullo_epi64(d_vec_0, alpha_vec);
            auto scaled_d_vec_1 = _mm512_mullo_epi64(d_vec_1, alpha_vec);
            auto scaled_d_vec_2 = _mm512_mullo_epi64(d_vec_2, alpha_vec);
            auto scaled_d_vec_3 = _mm512_mullo_epi64(d_vec_3, alpha_vec);
            auto scaled_d_vec_4 = _mm512_mullo_epi64(d_vec_4, alpha_vec);
            auto scaled_d_vec_5 = _mm512_mullo_epi64(d_vec_5, alpha_vec);
            auto scaled_d_vec_6 = _mm512_mullo_epi64(d_vec_6, alpha_vec);
            auto scaled_d_vec_7 = _mm512_mullo_epi64(d_vec_7, alpha_vec);

            // store
            uint64_t* store_loc = buf + row*buf_cols + col;
            _mm512_storeu_si512((void*)(store_loc), scaled_d_vec_0);
            _mm512_storeu_si512((void*)(store_loc + 8), scaled_d_vec_1);
            _mm512_storeu_si512((void*)(store_loc + 16), scaled_d_vec_2);
            _mm512_storeu_si512((void*)(store_loc + 24), scaled_d_vec_3);
            _mm512_storeu_si512((void*)(store_loc + 32), scaled_d_vec_4);
            _mm512_storeu_si512((void*)(store_loc + 40), scaled_d_vec_5);
            _mm512_storeu_si512((void*)(store_loc + 48), scaled_d_vec_6);
            _mm512_storeu_si512((void*)(store_loc + 56), scaled_d_vec_7);
        }
    }

    // first cleanup
    for (; col < col_num - 7; col+=8) {
        Packed_Loc col_loc = loc;

        // load packed vectors
        uint64_t* packed_loc_init = packed_data + col_loc.packed_num + col;
        auto p_vec_0 = _mm512_loadu_si512((void*)(packed_loc_init));

        for (int row = 0; row < row_num; row++) {
            // shift and mask
            auto sr_p_vec_0 = _mm512_srli_epi64(p_vec_0, col_loc.bit_offset);

            auto d_vec_0 = _mm512_and_si512(mask_vec, sr_p_vec_0);

            col_loc.bit_offset += num_bits;

            // check for overlap
            if (col_loc.bit_offset >= bits_type) {
                col_loc.packed_num += col_loc.eff_cols; 

                // load next vec
                uint64_t* packed_loc = packed_data + col_loc.packed_num + col;
                p_vec_0 = _mm512_loadu_si512((void*)(packed_loc));

                if (col_loc.bit_offset > bits_type) {
                    // shift and mask
                    int shift_amount = bits_type - (col_loc.bit_offset - num_bits);
                    auto sl_p_vec_0 = _mm512_slli_epi64(p_vec_0, shift_amount);

                    auto m_sl_p_vec_0 = _mm512_and_si512(mask_vec, sl_p_vec_0);

                    // add to result
                    d_vec_0 = _mm512_or_si512(d_vec_0, m_sl_p_vec_0);
                }

                col_loc.bit_offset &= (bits_type - 1);
            }

            auto scaled_d_vec_0 = _mm512_mullo_epi64(d_vec_0, alpha_vec);

            // store
            uint64_t* store_loc = buf + row*buf_cols + col;
            _mm512_storeu_si512((void*)(store_loc), scaled_d_vec_0);
        }
    }

    // final cleanup
    for (int row = 0; row < row_num; row++) {
        for (int col1 = col; col1 < col_num; col1++) {
            // get element
            buf[row*buf_cols + col1] = mask & (uint64_t)(packed_data[loc.packed_num + col1] >> loc.bit_offset);
        }

        // if overlaps to next
        loc.bit_offset += num_bits;

        if (loc.bit_offset >= bits_type) {
            loc.packed_num += loc.eff_cols;
            if (loc.bit_offset > bits_type) {
                for (int col1 = col; col1 < col_num; col1++) {
                    buf[row*buf_cols + col1] |= mask & (uint64_t)(packed_data[loc.packed_num + col1] << (bits_type - (loc.bit_offset - num_bits)));
                }
            }
            loc.bit_offset &= (bits_type - 1);
        }

        // scale
        for (int col1 = col; col1 < col_num; col1++) {
            buf[row*buf_cols + col1] *= alpha;
        }
    }
}
