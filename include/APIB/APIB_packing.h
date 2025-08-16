/**
 * @file APIB_unpacking.h
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains declaration of templated pack function
 * 
 */

#pragma once

#include "APIB/APIB_helper.h"

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
void pack(buf_type* buf, base_type* packed_data, Packed_Loc loc, int row_num, int col_num, int buf_cols, int num_bits);
