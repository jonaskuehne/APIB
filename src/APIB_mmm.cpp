/**
 * @file APIB_mmm.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains implementations of specializations of mini_mmm
 * 
 */

#include "APIB/APIB_mmm.h"
#include "APIB/APIB_helper.h"
#include <cstdint>
#include <immintrin.h>

// asm functions mini-mmm 8->32
extern "C" void mini_mmm_asm_8_32(uint8_t* a_buf, uint8_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U_4);
extern "C" void mini_mmm_asm_8_32_unr(uint8_t* a_buf, uint8_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U_4);

// asm functions mini-mmm 16->32
extern "C" void mini_mmm_asm_16_32(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_8(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_8_unr_32(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_4(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_4_unr_32(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_4_unr_64(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_2(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_2_unr_32(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_2_unr_64(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_2_unr_128(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_1(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_1_unr_32(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_1_unr_64(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_16_32_clean_1_unr_128(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, int k_inc, int P_U, int M_U);

// asm functions mini-mmm 64->64 using ifma
extern "C" void mini_mmm_asm_64_64_ifma(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_unr(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_4(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_4_unr_16(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_4_unr_32(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_2(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_2_unr_16(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_2_unr_32(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_2_unr_64(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_1(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_1_unr_16(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_1_unr_32(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_ifma_clean_1_unr_64(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);

// asm functions mini-mmm 64->64
extern "C" void mini_mmm_asm_64_64(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_4(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_4_unr_16(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_4_unr_32(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_2(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_2_unr_16(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_2_unr_32(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_1(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_1_unr_16(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_1_unr_32(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);
extern "C" void mini_mmm_asm_64_64_clean_1_unr_64(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int k_inc, int P_U, int M_U);

/**
 * @brief Computes mini-mmm of form C' <- C' + A'xB'
 * 
 * @tparam op_buf_type Type of buffer where A' and B' are unpacked to
 * @tparam res_buf_type Type of buffer where C' is unpacked to
 * @tparam version 0 if we can use IFMA for 64bit source, 1 if not
 * @param a_buf Pointer to buffer containing A
 * @param b_buf Pointer to buffer containing B
 * @param c_buf Pointer to buffer containing C
 * @param i_inc Increase in n dimension
 * @param j_inc Increase in m dimension
 * @param k_inc Increase in p dimension
 * @param M_U Block size in m dimension
 * @param P_U Block size in p dimension
 */
template <typename op_buf_type, typename res_buf_type, int version>
void mini_mmm(op_buf_type* a_buf, op_buf_type* b_buf, res_buf_type* c_buf, 
        int i_inc, int j_inc, int k_inc, int M_U, int P_U) {
    if (CHECK_MATRIX) printf("generic mini mmm\n");
    for (int i1 = 0; i1 < i_inc; i1++) {
        for (int j1 = 0; j1 < j_inc; j1++) {
            res_buf_type tmp = c_buf[i1*M_U + j1];
            for (int k1 = 0; k1 < k_inc; k1++) {
                tmp += (res_buf_type)(a_buf[i1*P_U + k1] * b_buf[k1*M_U + j1]); 
            }
            c_buf[i1*M_U + j1] = tmp;
        }
    }
}

// specialization for 16-32
template <>
void mini_mmm<uint16_t, uint32_t, 0>(uint16_t* a_buf, uint16_t* b_buf, uint32_t* c_buf, 
        int i_inc, int j_inc, int k_inc, int M_U, int P_U) {
    if (CHECK_MATRIX) printf("mini mmm with 16bit bufs and 32bit res\n");
    int i1 = 0;
    
    // main loop
    for(; i1 < i_inc - 15; i1+=16) {
        for (int j1 = 0; j1 < j_inc; j1+=16) {
            mini_mmm_asm_16_32(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
    }
    // everything cleanup from here

    for(; i1 < i_inc - 7; i1+=8) {
        int j1 = 0;
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_16_32_clean_8_unr_32(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=16) {
            mini_mmm_asm_16_32_clean_8(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
    }

    for(; i1 < i_inc - 3; i1+=4) {
        int j1 = 0;
        for (; j1 < j_inc - 63; j1+=64) {
            mini_mmm_asm_16_32_clean_4_unr_64(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_16_32_clean_4_unr_32(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc; j1+=16) {
            mini_mmm_asm_16_32_clean_4(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
    }

    for(; i1 < i_inc - 1; i1+=2) {
        int j1 = 0;
        for (; j1 < j_inc - 127; j1+=128) {
            mini_mmm_asm_16_32_clean_2_unr_128(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc - 63; j1+=64) {
            mini_mmm_asm_16_32_clean_2_unr_64(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_16_32_clean_2_unr_32(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc; j1+=16) {
            mini_mmm_asm_16_32_clean_2(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
    }

    for(; i1 < i_inc; i1+=1) {
        int j1 = 0;
        for (; j1 < j_inc - 127; j1+=128) {
            mini_mmm_asm_16_32_clean_1_unr_128(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc - 63; j1+=64) {
            mini_mmm_asm_16_32_clean_1_unr_64(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_16_32_clean_1_unr_32(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
        for (; j1 < j_inc; j1+=16) {
            mini_mmm_asm_16_32_clean_1(a_buf + i1*P_U, b_buf + (j1 << 1), c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        } 
    }
}

// specialization for 8-32
// buf[(row >> 2)*buf_cols + (col1 << 2) + (row & 3)]
template <>
void mini_mmm<uint8_t, uint32_t, 0>(uint8_t* a_buf, uint8_t* b_buf, uint32_t* c_buf, int i_inc, int j_inc, int k_inc, int M_U, int P_U) {
    if (CHECK_MATRIX) printf("mini mmm with 8bit bufs and 32bit res\n");
    int i1 = 0;
    int j1 = 0;

    // main loop
    for (; i1 < i_inc - 15; i1+=16) {
        for (j1=0; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_8_32_unr(a_buf + i1*P_U, b_buf + (j1 << 2), c_buf + i1*M_U + j1, k_inc, P_U, M_U << 2); 
        }
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_8_32(a_buf + i1*P_U, b_buf + (j1 << 2), c_buf + i1*M_U + j1, k_inc, P_U, M_U << 2); 
        }
    }

    // cleanup: i_inc not div by 16
    if (i1 < i_inc) {
        // change tile config
        __tilecfg tile_data;
        _tile_storeconfig(&tile_data);
        int rows = i_inc - i1;

        tile_data.rows[0] = rows;
        tile_data.rows[1] = rows;
        tile_data.rows[3] = rows;
        _tile_loadconfig (&tile_data);

        for (j1=0; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_8_32_unr(a_buf + i1*P_U, b_buf + (j1 << 2), c_buf + i1*M_U + j1, k_inc, P_U, M_U << 2); 
        }

        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_8_32(a_buf + i1*P_U, b_buf + (j1 << 2), c_buf + i1*M_U + j1, k_inc, P_U, M_U << 2); 
        }

        // restore tile config
        tile_data.rows[0] = MAX_ROWS;
        tile_data.rows[1] = MAX_ROWS;
        tile_data.rows[3] = MAX_ROWS;
        _tile_loadconfig (&tile_data);
    }

    // cleanup: j_inc not div by 16
    if (j1 < j_inc) {
        // change tile config
        __tilecfg tile_data;
        _tile_storeconfig(&tile_data);
        int colsb = (j_inc - j1) << 2;

        tile_data.colsb[0] = colsb;
        tile_data.colsb[1] = colsb;
        tile_data.colsb[4] = colsb;
        tile_data.colsb[5] = colsb;
        _tile_loadconfig (&tile_data);

        for (i1=0; i1 < i_inc - 15; i1+=16) {
            mini_mmm_asm_8_32(a_buf + i1*P_U, b_buf + (j1 << 2), c_buf + i1*M_U + j1, k_inc, P_U, M_U << 2); 
        }

        // restore tile config
        tile_data.colsb[0] = MAX_COLSB;
        tile_data.colsb[1] = MAX_COLSB;
        tile_data.colsb[4] = MAX_COLSB;
        tile_data.colsb[5] = MAX_COLSB;
        _tile_loadconfig (&tile_data);
    }

    // cleanup: i_inc and j_inc not div by 16
    if (i1 < i_inc && j1 < j_inc) {
        // change tile config
        __tilecfg tile_data;
        _tile_storeconfig(&tile_data);
        int rows = i_inc - i1;
        int colsb = (j_inc - j1) << 2;

        tile_data.rows[0] = rows;
        tile_data.rows[1] = rows;
        tile_data.rows[3] = rows;

        tile_data.colsb[0] = colsb;
        tile_data.colsb[1] = colsb;
        tile_data.colsb[4] = colsb;
        tile_data.colsb[5] = colsb;

        _tile_loadconfig (&tile_data);

        mini_mmm_asm_8_32(a_buf + i1*P_U, b_buf + (j1 << 2), c_buf + i1*M_U + j1, k_inc, P_U, M_U << 2); 

        // restore tile config
        tile_data.rows[0] = MAX_ROWS;
        tile_data.rows[1] = MAX_ROWS;
        tile_data.rows[3] = MAX_ROWS;

        tile_data.colsb[0] = MAX_COLSB;
        tile_data.colsb[1] = MAX_COLSB;
        tile_data.colsb[4] = MAX_COLSB;
        tile_data.colsb[5] = MAX_COLSB;

        _tile_loadconfig (&tile_data);
    }
}

// specialization for 64->64 with ifma
template <>
void mini_mmm<uint64_t, uint64_t, 0>(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int i_inc, int j_inc, int k_inc, int M_U, int P_U) {
    if (CHECK_MATRIX) printf("64bit mini mmm with ifma\n");
    int i1 = 0;

    // main loop
    for (; i1 < i_inc - 7; i1+=8) {
        int j1 = 0; 
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_64_64_ifma_unr(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64_ifma(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }
    // everything cleanup from here

    for (; i1 < i_inc - 3; i1+=4) {
        int j1 = 0;
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_64_64_ifma_clean_4_unr_32(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_64_64_ifma_clean_4_unr_16(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64_ifma_clean_4(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }

    for (; i1 < i_inc - 1; i1+=2) {
        int j1 = 0;    
        for (; j1 < j_inc - 63; j1+=64) {
            mini_mmm_asm_64_64_ifma_clean_2_unr_64(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_64_64_ifma_clean_2_unr_32(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_64_64_ifma_clean_2_unr_16(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64_ifma_clean_2(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }

    for (; i1 < i_inc; i1+=1) {
        int j1 = 0; 
        for (; j1 < j_inc - 63; j1+=64) {
            mini_mmm_asm_64_64_ifma_clean_1_unr_64(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_64_64_ifma_clean_1_unr_32(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_64_64_ifma_clean_1_unr_16(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64_ifma_clean_1(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }
}

template <>
void mini_mmm<uint64_t, uint64_t, 1>(uint64_t* a_buf, uint64_t* b_buf, uint64_t* c_buf, int i_inc, int j_inc, int k_inc, int M_U, int P_U) {
    if (CHECK_MATRIX) printf("64bit mini mmm\n");
    int i1 = 0;

    // main loop
    for (; i1 < i_inc - 7; i1+=8) {
        int j1 = 0; 
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }
    // everything cleanup from here

    for (; i1 < i_inc - 3; i1+=4) {
        int j1 = 0;
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_64_64_clean_4_unr_32(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_64_64_clean_4_unr_16(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64_clean_4(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }
    for (; i1 < i_inc; i1+=2) {
        int j1 = 0;    
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_64_64_clean_2_unr_32(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_64_64_clean_2_unr_16(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64_clean_2(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }

    for (; i1 < i_inc; i1+=1) {
        int j1 = 0; 
        for (; j1 < j_inc - 63; j1+=64) {
            mini_mmm_asm_64_64_clean_1_unr_64(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 31; j1+=32) {
            mini_mmm_asm_64_64_clean_1_unr_32(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc - 15; j1+=16) {
            mini_mmm_asm_64_64_clean_1_unr_16(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
        for (; j1 < j_inc; j1+=8) {
            mini_mmm_asm_64_64_clean_1(a_buf + i1*P_U, b_buf + j1, c_buf + i1*M_U + j1, k_inc, P_U, M_U);
        }
    }
}
