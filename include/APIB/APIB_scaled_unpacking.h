/**
 * @file APIB_scaled_unpacking.h
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains declarations for templated scaled unpack functions
 * 
 */

#pragma once

#include "APIB/APIB_helper.h"

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
void scaled_unpack(buf_type* buf, base_type* packed_data, buf_type alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits);

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
void scaled_amx_unpack(buf_type* buf, base_type* packed_data, buf_type alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits);

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
template <typename buf_type, typename base_type>
void scaled_madd_unpack(buf_type* buf, base_type* packed_data, buf_type alpha, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits);
