/**
 * @file APIB_mmm.h
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains declarations for templated mini-mmm function
 * 
 */

#pragma once

#include <cstdint>

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
template<typename op_buf_type, typename res_buf_type, int version>
void mini_mmm(op_buf_type* a_buf, op_buf_type* b_buf, res_buf_type* c_buf, 
        int i_inc, int j_inc, int k_inc, int M_U, int P_U);
