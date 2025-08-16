/**
 * @file APIB.h
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains struct definitions for APIB matrices and their multiplication
 * 
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <variant>
#include <vector>
#include "APIB/APIB_helper.h"
#include "APIB/APIB_packing.h"
#include "APIB/APIB_unpacking.h"
#include "APIB/APIB_scaled_unpacking.h"
#include "APIB/APIB_mmm.h"
#include "APIB/APIB_consts.h"

/**
 * @brief Generic matrix object for APIB
 * 
 * @tparam num_bits Number of bits per element
 * @tparam num_rows Number of rows of matrix
 * @tparam num_cols Number of columns of matrix
 * @tparam input_type Type of buffer with initialization data
 */
template <int num_bits, int num_rows, int num_cols, typename input_type>
struct APIB_Mat {
    // type in which the packed elements are stored
    using base_type = typename get_base_type<num_bits>::type;

    const static int bits_type = 8 * sizeof(base_type);
    // will need this later, const for efficiency
    const static int bits_type_exp = Ceil_Log2<bits_type>::value;
    // how many base_type elements needed per column
    const static int num_packed_els_per_col = 1 + (num_bits*num_rows + (bits_type - 1)) / bits_type; 
    // total amount of element used in packed format
    const static int num_packed_els = num_cols * num_packed_els_per_col;
    // blockwidth, benchmarked parameter
    // need exp for efficient computation later on
    const int block_width_exp = APIB_BLOCK_EXP;
    const int block_width = 1 << block_width_exp;
    // for slices
    static constexpr int num_bits_value = num_bits;
    static constexpr int num_rows_value = num_rows;

    /**
     * @brief Computes the size matrix in the packed format in bytes
     * 
     * @return int Size in bytes
     */
    static int get_size() {
        return num_packed_els * sizeof(base_type);
    }

    /**
     * @brief Construct a new uninitialized APIB Mat object, internal data points to NULL.
     * To initialize (create and zero data) call reset() on the object
     * 
     */
    APIB_Mat() : packed_data(NULL) {
        // check size
        if (num_bits > 64) {
            printf("APIB: trying to create APIB_Mat object with %d bits, but 64 bit is max\n", num_bits);
            if (CRASH_IF_ERROR) exit(1);
        }

        // check buffer type
        if (num_bits > 8 * sizeof(input_type)) {
            printf("APIB: trying to init APIB_Mat object with %d bits, but buffer can hold at most %ld-bit elements\n", 
                    num_bits, 8 * sizeof(input_type));
            if (CRASH_IF_ERROR) exit(1);
        }
    };

    /**
     * @brief Construct a initialized new APIB Mat object
     * 
     * @param input_data Input data, expected to be non-null and in row-major format
     */
    APIB_Mat(input_type* input_data) : APIB_Mat() {
        // initialized buffer
        reset();
        // initialize block-wise
        #pragma omp parallel for shared(input_data, packed_data) default(none) schedule(static)
        for (int col_block = 0; col_block < num_cols; col_block+=block_width) {
            int col_lim = std::min(num_cols - col_block, block_width);

            // use packing utility
            pack<input_type, base_type>(input_data + col_block, packed_data, get_loc(0, col_block), num_rows, col_lim, num_cols, num_bits); 
        }
    }

    /**
     * @brief Unpacks the packed data to a buffer
     * 
     * @tparam buf_type Type of the destination buffer, expected to fit the elements
     * @param buf Pointer to buffer
     * @param buf_cols Number of columns of the buffer
     */
    template <typename buf_type>
    void unpack_data(buf_type* buf, int buf_cols) {
        // check size
        if (num_bits > 8 * sizeof(buf_type)) {
            printf("APIB: trying to unpack APIB_Mat object with %d bits, but buffer can hold at most %ld-bit elements\n", 
                    num_bits, 8 * sizeof(buf_type));
            if (CRASH_IF_ERROR) exit(1);
        }

        // unpack block-wise
        #pragma omp parallel for shared(buf, packed_data, buf_cols) default(none) schedule(static)
        for (int col_block = 0; col_block < num_cols; col_block+=block_width) {
            int col_lim = std::min(num_cols - col_block, block_width);

            // use unpacking utility
            unpack<buf_type, base_type>(buf + col_block, packed_data, get_loc(0, col_block), num_rows, col_lim, buf_cols, num_bits); 
        }
    }

    /**
     * @brief Resets the buffer holding the packed data or creates it if necessary
     * 
     */
    void reset() {
        if (packed_data != NULL) {
            free(packed_data);
        }

        packed_data = (base_type*)aligned_alloc(4096, num_packed_els*sizeof(base_type));
        // need to zero memory for the packing process
        // no aligned calloc, hence need to do in two steps
        memset((void*)packed_data, 0, sizeof(base_type) * num_packed_els);
    }

    /**
     * @brief Copy Constructor (deep copy), if other.packed_data is null, so is the one of the resulting matrix
     * 
     * @param other Source Object
     */
    APIB_Mat(const APIB_Mat& other) : packed_data((base_type*)aligned_alloc(4096, other.num_packed_els*sizeof(base_type))) {
        if (packed_data != NULL) {
            free(packed_data);
        }

        if (other.packed_data != NULL) {
            // initialize to size of source object and copy data
            packed_data = (base_type*)aligned_alloc(4096, other.num_packed_els*sizeof(base_type));
            memcpy(packed_data, other.packed_data, other.num_packed_els * sizeof(base_type));
        } else {
            packed_data = NULL;
        }
    }

    /**
     * @brief Destroy the APIB Mat, frees all memory
     * 
     */
    ~APIB_Mat() {
        if (packed_data != NULL) {
            free(packed_data);
        }
    }

    /**
     * @brief Copy Assignment (deep copy), if other.packed_data is null, so is the one of the resulting matrix. 
     * 
     * @param other Source Object
     * @return APIB_Mat& 
     */
    APIB_Mat& operator=(const APIB_Mat& other) {
        if (this == &other) {
            return *this;
        }

        if (packed_data != NULL) {
            free(packed_data);
        }

        if (other.packed_data != NULL) {
            // initialize to size of source object and copy data
            packed_data = (base_type*)aligned_alloc(4096, other.num_packed_els*sizeof(base_type));
            memcpy(packed_data, other.packed_data, other.num_packed_els * sizeof(base_type));
        } else {
            packed_data = NULL;
        }

        return *this;
    }

    /**
     * @brief Prints all elements of the matrix object
     * 
     */
    void print_els() const {
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                printf("%lu\t", get_el<uint64_t>(i, j));
            }
            printf("\n");
        }
    }

    /**
     * @brief Creates a location of the element at (row, column)
     * 
     * @param row Row
     * @param col Column
     * @return Packed_Loc 
     */
    Packed_Loc get_loc(int row, int col) const {
        // block number
        int block = col >> block_width_exp;
        // offset of block (in packed data)
        int block_offset = (block << block_width_exp) * num_packed_els_per_col;
        // effective columns in block
        int eff_cols = std::min(num_cols - (block << block_width_exp), block_width);

        // bit number in column
        int bit_num = row * num_bits;
        // offset in packed element
        int bit_offset = bit_num & (bits_type - 1);
        // row in packed format
        int packed_row = bit_num >> bits_type_exp;
        // position in packed format
        int packed_num = block_offset + packed_row*eff_cols + (col & (block_width - 1));

        // create element
        Packed_Loc loc(packed_num, eff_cols, bit_offset);
        return loc;
    }

    /**
     * @brief Returns a single element at pos (row, col), meant for debugging purposes
     * 
     * @tparam t Type of returned element
     * @param row Row
     * @param col Column
     * @return t 
     */
    template<typename t>
    t get_el(int row, int col) const {
        // get location
        Packed_Loc loc = get_loc(row, col); 

        // set mask
        base_type mask;
        if (num_bits >= 64) {
            mask = (base_type)(-1);
        } else if (num_bits >= 63) {
            mask = (base_type)(~(1L << num_bits));
        } else {
            mask = (base_type)((1L << num_bits) - 1);
        }

        // normal case 
        int64_t res = mask & (packed_data[loc.packed_num] >> loc.bit_offset);
        // overlap to next
        if (loc.bit_offset + num_bits > bits_type) {
            int overlap = num_bits - (bits_type - loc.bit_offset);
            res |= mask & (packed_data[loc.packed_num + loc.eff_cols] << (num_bits - overlap));
        }

        return (t)res;
    }

    // pointer to internal data representation
    base_type* packed_data;
};

/**
 * @brief Matrix object for APIB containing of horizontal slices, in a slice all 
 * elements take the same number of bits to store them. The slices are represented
 * as a vector of APIB_Mat objects.
 * 
 * @tparam cols Number of columns of the slices, equal for all
 * @tparam Specs List of APIB_Slice_Specs, describing each slice
 */
template <int cols, typename... Specs>
struct APIB_Slice_Mat {
    template <typename Spec>
    // list of types of the slices
    using Mat_Type = APIB_Mat<Spec::bits, Spec::rows, cols, typename Spec::Type>;
    // removed duplicates, variant can't handle duplicates
    using UniqueTypes = typename unique_types<Mat_Type<Specs>...>::type;
    // convert to variant
    using VariantType = typename tuple_to_variant<UniqueTypes>::type;
    // vector containing the slices
    std::vector<VariantType> slices;
    /**
     * @brief Construct an initialized APIB_Slice_Mat object
     * 
     * @param input_buffers Buffers to initialized the slices
     */
    APIB_Slice_Mat(typename Specs::Type*... input_buffers) : slices({ Mat_Type<Specs>(input_buffers)... }){}
    /**
     * @brief Construct an APIB_Slice_Mat object without any slices
     * 
     */
    APIB_Slice_Mat() : slices({}) {}

    /**
     * @brief Copy-Constructor (deep copy)
     * 
     * @param other Source Object
     */
    APIB_Slice_Mat(const APIB_Slice_Mat& other) : slices(other.slices) {}

    /**
     * @brief Copy-Assignment (deep copy)
     * 
     * @param other Source Object
     * @return APIB_Slice_Mat&
     */
    APIB_Slice_Mat& operator=(const APIB_Slice_Mat& other) {
        // nothing has to be done for self assignment
        if (this == &other) {
            return *this;
        }

        slices = other.slices;

        return *this;
    }

    /**
     * @brief Computes the size matrix in the packed format in bytes
     * 
     * @return int Size in bytes
     */
    int get_size() {
        // compute sum over all slices
        int size = 0;
        int num_slices = slices.size();
        for (int i = 0; i < num_slices; i++) {
            // visit to access slice
            std::visit([&size](auto& a) {
                size += a.get_size();
            }, slices[i]); 
        }
        return size;
    }
};

/**
 * @brief Utility function to compute a schedule with which matrix slices a buffer
 * should be filled. Intended for internal use
 * 
 * @tparam N_U Block size
 * @tparam n Number of rows in whole A matrix
 * @tparam p Number of columns of the slices
 * @tparam Specs List of APIB_Slice_Specs describing the slices
 * @param A APIB_Slice_Matrix Object
 * @param index_vec Vector to hold slice indices
 * @param row_vec Vector to hold number of rows to unpack
 * @param mat_offset_vec Vector to hold offset of rows in slice
 * @param buffer_offset_vec Vector to hold offset of rows in buffer
 */
template<int N_U, int n, int p, typename... Specs>
void build_slice_vecs(const APIB_Slice_Mat<p, Specs...>& A, 
    std::vector<std::vector<int>>& index_vec, 
    std::vector<std::vector<int>>& row_vec, 
    std::vector<std::vector<int>>& mat_offset_vec, 
    std::vector<std::vector<int>>& buffer_offset_vec) {

    // build vecs describing how to load matrices
    int seq_ind = 0;
    int rest = 0;
    int mat_offset = 0;
    int buffer_offset = 0;
    int i_inc = N_U;
    index_vec = {};
    row_vec = {};
    mat_offset_vec = {};
    buffer_offset_vec = {};
    for(int i = 0; i < n; i+=i_inc) {
	buffer_offset = 0;
        int i_inc = std::min(N_U, n - i);
        std::vector<int> indices = {};
        std::vector<int> rows = {};
        std::vector<int> mat_offsets = {};
        std::vector<int> buffer_offsets = {};

        // rest of last iteration enough for this one
        if (rest > i_inc) {
            indices.push_back(seq_ind);
            rows.push_back(i_inc);
            mat_offsets.push_back(mat_offset);
            buffer_offsets.push_back(buffer_offset);

            rest -= i_inc;
            mat_offset += i_inc;
            buffer_offset += i_inc;
        // need atleast another matrix
        } else {
            int curr = 0;
            // have something left from last iteration
            if (rest != 0) {
                indices.push_back(seq_ind);
                rows.push_back(rest);
                mat_offsets.push_back(mat_offset);
                buffer_offsets.push_back(buffer_offset);
                
                seq_ind++;
                buffer_offset += rest;
                curr = rest;
                mat_offset = 0;
                rest = 0;
            }
            
            // fill
            while (curr < i_inc) {
                // use visit and functor to access slice
                std::visit([&seq_ind, &rest, &curr, i_inc, &mat_offset, 
                        &buffer_offset, &indices, &rows, &mat_offsets, 
                        &buffer_offsets](const auto& slice) {
                    const int num_rows = std::decay_t<decltype(slice)>::num_rows_value;
                    curr += num_rows;
                    if (curr > i_inc) {
                        rest = curr - i_inc; 
                        int amnt = num_rows - rest;
                        indices.push_back(seq_ind);
                        rows.push_back(amnt);
                        mat_offsets.push_back(mat_offset);
                        buffer_offsets.push_back(buffer_offset);

                        buffer_offset += amnt;
                        mat_offset = amnt;
                    } else {
                        indices.push_back(seq_ind);
                        rows.push_back(num_rows);
                        mat_offsets.push_back(mat_offset);
                        buffer_offsets.push_back(buffer_offset);

                        buffer_offset += num_rows;
                        seq_ind++;
                    }
                }, A.slices[seq_ind]);
            }
        }
        index_vec.push_back(indices);
        row_vec.push_back(rows);
        mat_offset_vec.push_back(mat_offsets);
        buffer_offset_vec.push_back(buffer_offsets);
    }
}

/**
 * @brief Performs matrix multiplication D <- a*C + b*AxB where B, C, D are normal
 * APIB matrices and A is a sliced APIB matrix
 * 
 * @tparam num_bits_alpha Number of bits of a
 * @tparam num_bits_beta Number of bits of b
 * @param A APIB Slice matrix A
 * @param B APIB matrix B
 * @param C APIB matrix C
 * @param alpha Scaling factor of C
 * @param beta Scaling factor of AxB
 * @return APIB matrix D
 */
template <int num_bits_alpha, int num_bits_beta, 
         int num_bits_B, int num_bits_C, int n, int m, int p, 
         typename type_B, typename type_C, typename... Specs>
auto mmm(const APIB_Slice_Mat<p, Specs...>& A, 
        const APIB_Mat<num_bits_B, p, m, type_B>& B, 
        const APIB_Mat<num_bits_C, n, m, type_C>& C,
        const uint64_t alpha, const uint64_t beta) {
    
    // make sure A has the same number of rows as C (and same number of cols as B)
    const int sum_rows = get_total_rows<Specs...>::value;
    static_assert(sum_rows == n);
    // D has to fit elements from both parts of A
    const int num_bits_A = get_max_bits<Specs...>::value;
    // number of bits result and operands will be (with scaling)
    const int res_bits = Min_Bits_Res<num_bits_alpha, num_bits_beta, num_bits_A, num_bits_B, num_bits_C, p>::value;
    const int op_bits = Min_Bits_Op<num_bits_beta, num_bits_A, num_bits_B>::value;

    using buf_type = typename get_op_buf_type<op_bits, res_bits>::type;
    using res_buf_type = typename get_res_buf_type<op_bits, res_bits>::type;


    // initialize matrix object
    APIB_Mat<res_bits, n, m, res_buf_type> D;  
    D.reset();

    using base_type_B = typename std::decay_t<decltype(B)>::base_type;
    using base_type_C = typename std::decay_t<decltype(C)>::base_type;
    using base_type_D = typename std::decay_t<decltype(D)>::base_type;

    // only relevant for 64bit results, if result at most 52 bits we can use ifma instructions
    const int mmm_version = num_bits_A + num_bits_B + num_bits_beta > 52;

    // get blocksizes
    const int N_U = Get_N_U<8*sizeof(buf_type), mmm_version>::value;
    const int M_U = Get_M_U<8*sizeof(buf_type), mmm_version>::value;
    const int P_U = Get_P_U<8*sizeof(buf_type), mmm_version>::value;

    // prefer to scale A, more efficient as unpacking of B sometimes more complex
    const int scale_A = num_bits_A + num_bits_beta <= 8 * sizeof(buf_type); 

    // check if we use special instructions
    const int use_amx = sizeof(buf_type) == 1 && sizeof(res_buf_type) == 4;
    const int use_madd = sizeof(buf_type) == 2 && sizeof(res_buf_type) == 4;

    // debug output
    if (CHECK_RES_BITS) printf("needs %d bits\n", res_bits);
    if (CHECK_TYPES) printf("Operands: %lu\nResult: %lu, Version:%d\n", 
            8*sizeof(buf_type), 8*sizeof(res_buf_type), mmm_version);

    // vectors to hold information on how to process slices
    std::vector<std::vector<int>> index_vec;
    std::vector<std::vector<int>> row_vec;
    std::vector<std::vector<int>> mat_offset_vec;
    std::vector<std::vector<int>> buffer_offset_vec;
    // fill
    build_slice_vecs<N_U, n>(A, index_vec, row_vec, mat_offset_vec, buffer_offset_vec);

    #pragma omp parallel shared(A, B, C, D, alpha, beta, index_vec, row_vec, mat_offset_vec, buffer_offset_vec) default(none)
    {

        // setup amx, needs to be done for all threads
        if (use_amx) {
            // get permission to use AMX
            if (!set_tiledata_use()) exit(-1);
            __tilecfg tile_data = {0};
            // load tile configuration 
            init_tile_config (&tile_data, MAX_ROWS, MAX_COLSB);
        }

        // allocate buffers
        buf_type* a_buf = (buf_type*) aligned_alloc(4096, N_U*P_U*sizeof(buf_type)); 
        buf_type* b_buf = (buf_type*) aligned_alloc(4096, P_U*M_U*sizeof(buf_type)); 
        res_buf_type* c_buf = (res_buf_type*) aligned_alloc(4096, N_U*M_U*sizeof(res_buf_type)); 

        // initialize inc steps for all dimensions
        int i_inc = N_U;
        int j_inc = M_U;
        int k_inc = P_U;


        // MMM
        // swap i and j loop to make sure no race condition when writing to D
        #pragma omp for schedule(static)
        for(int j = 0; j < m; j+=j_inc) {
            // make sure no out of bounds memory accesses
            j_inc = (m-j > M_U) ? M_U : m - j;

            int vec_index = 0;
            for(int i = 0; i < n; i+=i_inc, vec_index++) {
                // make sure no out of bounds memory accesses
                i_inc = (n-i > N_U) ? N_U : n - i;

                // unpack and scale C into continuous memory
                scaled_unpack<res_buf_type, base_type_C>
                        (c_buf, C.packed_data, (res_buf_type)alpha, 
                        C.get_loc(i, j), i_inc, j_inc, M_U, num_bits_C);

                for(int k = 0; k < p; k+=k_inc) {
                    // make sure no out of bounds memory accesses
                    k_inc = (p-k > P_U) ? P_U : p - k;
                    
                    // scale A
                    if (scale_A) {
                        // unpack and scale A into continuous memory
                        int vec_length = index_vec[vec_index].size();
                        // go through vector describing what and how to unpack
                        for (int vec_entry = 0; vec_entry < vec_length; vec_entry++) {
                            // visit to access slice
                            std::visit([a_buf, beta, i, k, i_inc, k_inc, P_U, 
                                    &mat_offset_vec, &row_vec, &buffer_offset_vec, 
                                    vec_index, vec_entry](auto& a) {
                                // extract base type (32 or 64 bit) of slice
                                using base_type_slice = typename std::decay_t<decltype(a)>::base_type;
                                const int num_bits_slice = std::decay_t<decltype(a)>::num_bits_value;
                                // actual unpack and scale with parameters from vector
                                scaled_unpack<buf_type, base_type_slice>
                                    (a_buf + P_U * buffer_offset_vec[vec_index][vec_entry], a.packed_data, 
                                     (buf_type)beta, a.get_loc(mat_offset_vec[vec_index][vec_entry], k), 
                                    row_vec[vec_index][vec_entry], k_inc, P_U, num_bits_slice);
                            }, A.slices[index_vec[vec_index][vec_entry]]); 
                        }
                        // unpack B into continuous memory
                        // can use amx
                        if (use_amx) {
                            // check if block not full, need to zero some memory in this case
                            if (k_inc != P_U) {
                                int num_bytes_used = (k_inc >> 2) * (M_U << 2);
                                memset((void*)(b_buf + num_bytes_used), 0, 
                                        P_U*M_U - num_bytes_used);
                            }
                            // unpack to format for AMX
                            amx_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 2, num_bits_B);
                        // use madd (VNNI)
                        } else if (use_madd) {
                            madd_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 1, num_bits_B);
                        // normal format (64bit source)
                        } else {
                            unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U, num_bits_B);
                        }
                    // scale B
                    } else {
                        // unpack A into continuous memory
                        int vec_length = index_vec[vec_index].size();
                        // go through vector describing what and how to unpack
                        for (int vec_entry = 0; vec_entry < vec_length; vec_entry++) {
                            // visit to access slice
                            std::visit([a_buf, i, k, i_inc, k_inc, P_U, 
                                    &mat_offset_vec, &row_vec, &buffer_offset_vec, 
                                    vec_index, vec_entry](auto& a) {
                                // extract base type (32 or 64 bit) of slice
                                using base_type_slice = typename std::decay_t<decltype(a)>::base_type;
                                const int num_bits_slice = std::decay_t<decltype(a)>::num_bits_value;
                                // actual unpack with parameters from vector
                                unpack<buf_type, base_type_slice>
                                    (a_buf + P_U * buffer_offset_vec[vec_index][vec_entry], a.packed_data, 
                                     a.get_loc(mat_offset_vec[vec_index][vec_entry], k), 
                                    row_vec[vec_index][vec_entry], k_inc, P_U, num_bits_slice);
                            }, A.slices[index_vec[vec_index][vec_entry]]); 
                        }
                        // unpack and scale B into continuous memory
                        // can use amx
                        if (use_amx) {
                            // check if block not full, need to zero some memory in this case
                            if (k_inc != P_U) {
                                int num_bytes_used = (k_inc >> 2) * (M_U << 2);
                                memset((void*)(b_buf + num_bytes_used), 0, P_U*M_U - num_bytes_used);
                            }
                            // scale and unpack to AMX-format
                            scaled_amx_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, (buf_type)beta, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 2, num_bits_B);
                        // use m_add (VNNI)
                        } else if (use_madd) {
                            scaled_madd_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, (buf_type)beta, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 1, num_bits_B);
                        // normal format (64bit source)
                        } else {
                            scaled_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, (buf_type)beta, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U, num_bits_B);
                        }
                    }
                    
                    // execute mini-mmm on buffers with unpacked data
                    mini_mmm<buf_type, res_buf_type, mmm_version>(a_buf, b_buf, c_buf, 
                            i_inc, j_inc, k_inc, M_U, P_U);
                }
                
                // pack result to D
                pack<res_buf_type, base_type_D>
                        (c_buf, D.packed_data, D.get_loc(i, j), 
                        i_inc, j_inc, M_U, res_bits);
            }
        }
        
        // used buffers
        free(a_buf);
        free(b_buf);
        free(c_buf);

        // reset tiles after use
        if (use_amx) _tile_release();
    }
    return D;
}

/**
 * @brief Performs matrix multiplication D <- a*C + b*AxB where A, B, C, D are normal
 * APIB matrices
 * 
 * @tparam num_bits_alpha Number of bits of a
 * @tparam num_bits_beta Number of bits of b
 * @param A APIB matrix A
 * @param B APIB matrix B
 * @param C APIB matrix C
 * @param alpha Scaling factor of C
 * @param beta Scaling factor of AxB
 * @return APIB matrix D
 */
template <int num_bits_alpha, int num_bits_beta, 
         int num_bits_A, int num_bits_B, int num_bits_C, 
         int n, int m, int p, 
         typename type_A, typename type_B, typename type_C>
auto mmm(const APIB_Mat<num_bits_A, n, p, type_A>& A, 
        const APIB_Mat<num_bits_B, p, m, type_B>& B, 
        const APIB_Mat<num_bits_C, n, m, type_C>& C,
        const uint64_t alpha, const uint64_t beta) {

    // number of bits result and operands will be (with scaling)
    const int res_bits = Min_Bits_Res<num_bits_alpha, num_bits_beta, num_bits_A, num_bits_B, num_bits_C, p>::value;
    const int op_bits = Min_Bits_Op<num_bits_beta, num_bits_A, num_bits_B>::value;

    using buf_type = typename get_op_buf_type<op_bits, res_bits>::type;
    using res_buf_type = typename get_res_buf_type<op_bits, res_bits>::type;

    // initialize matrix object
    APIB_Mat<res_bits, n, m, res_buf_type> D;  
    D.reset();

    using base_type_A = typename std::decay_t<decltype(A)>::base_type;
    using base_type_B = typename std::decay_t<decltype(B)>::base_type;
    using base_type_C = typename std::decay_t<decltype(C)>::base_type;
    using base_type_D = typename std::decay_t<decltype(D)>::base_type;

    // only relevant for 64bit results, if result at most 52 bits we can use ifma instructions
    const int mmm_version = num_bits_A + num_bits_B + num_bits_beta > 52;

    // get blocksizes
    const int N_U = Get_N_U<8*sizeof(buf_type), mmm_version>::value;
    const int M_U = Get_M_U<8*sizeof(buf_type), mmm_version>::value;
    const int P_U = Get_P_U<8*sizeof(buf_type), mmm_version>::value;

    // prefer to scale A, more efficient as unpacking of B sometimes more complex
    const int scale_A = num_bits_A + num_bits_beta <= 8 * sizeof(buf_type); 

    // check if we use special instructions
    const int use_amx = sizeof(buf_type) == 1 && sizeof(res_buf_type) == 4;
    const int use_madd = sizeof(buf_type) == 2 && sizeof(res_buf_type) == 4;

    // debug output
    if (CHECK_RES_BITS) printf("needs %d bits\n", res_bits);
    if (CHECK_TYPES) printf("Operands: %lu\nResult: %lu\n", 8*sizeof(buf_type), 8*sizeof(res_buf_type));

    #pragma omp parallel shared(A, B, C, D, alpha, beta) default(none)
    {

        // setup amx, needs to be done for all threads
        if (use_amx) {
            // get permission to use AMX
            if (!set_tiledata_use()) exit(-1);
            __tilecfg tile_data = {0};
            // load tile configuration 
            init_tile_config (&tile_data, MAX_ROWS, MAX_COLSB);
        }

        // allocate buffers
        buf_type* a_buf = (buf_type*) aligned_alloc(4096, N_U*P_U*sizeof(buf_type)); 
        buf_type* b_buf = (buf_type*) aligned_alloc(4096, P_U*M_U*sizeof(buf_type)); 
        res_buf_type* c_buf = (res_buf_type*) aligned_alloc(4096, N_U*M_U*sizeof(res_buf_type)); 

        // initialize inc steps for all dimensions
        int i_inc = N_U;
        int j_inc = M_U;
        int k_inc = P_U;

        // MMM
        // swap i and j loop to make sure no race condition when writing to D
        #pragma omp for schedule(static)
        for(int j = 0; j < m; j+=j_inc) {
            // make sure no out of bounds memory accesses
            j_inc = (m-j > M_U) ? M_U : m - j;
            for(int i = 0; i < n; i+=i_inc) {
                // make sure no out of bounds memory accesses
                i_inc = (n-i > N_U) ? N_U : n - i;

                // unpack and scale C into continuous memory
                scaled_unpack<res_buf_type, base_type_C>
                        (c_buf, C.packed_data, (res_buf_type)alpha, 
                        C.get_loc(i, j), i_inc, j_inc, M_U, num_bits_C);

                for(int k = 0; k < p; k+=k_inc) {
                    // make sure no out of bounds memory accesses
                    k_inc = (p-k > P_U) ? P_U : p - k;
                    
                    // scale A
                    if (scale_A) {
                        // unpack and scale A into continuous memory
                        scaled_unpack<buf_type, base_type_A>
                                (a_buf, A.packed_data, (buf_type)beta, A.get_loc(i, k), 
                                i_inc, k_inc, P_U, num_bits_A);
                        // unpack B into continuous memory
                        // can use amx
                        if (use_amx) {
                            // check if block not full, need to zero some memory in this case
                            if (k_inc != P_U) {
                                int num_bytes_used = (k_inc >> 2) * (M_U << 2);
                                memset((void*)(b_buf + num_bytes_used), 0, 
                                        P_U*M_U - num_bytes_used);
                            }
                            // unpack to format for AMX
                            amx_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 2, num_bits_B);
                        // use madd (VNNI)
                        } else if (use_madd) {
                            madd_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 1, num_bits_B);
                        // normal format (64bit source)
                        } else {
                            unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U, num_bits_B);
                        }
                    // scale B
                    } else {
                        // unpack A into continuous memory
                        unpack<buf_type, base_type_A>
                                (a_buf, A.packed_data, A.get_loc(i, k), 
                                i_inc, k_inc, P_U, num_bits_A);
                        // unpack and scale B into continuous memory
                        // can use amx
                        if (use_amx) {
                            // check if block not full, need to zero some memory in this case
                            if (k_inc != P_U) {
                                int num_bytes_used = (k_inc >> 2) * (M_U << 2);
                                memset((void*)(b_buf + num_bytes_used), 0, P_U*M_U - num_bytes_used);
                            }
                            // scale and unpack to AMX-format
                            scaled_amx_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, (buf_type)beta, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 2, num_bits_B);
                        // use m_add (VNNI)
                        } else if (use_madd) {
                            scaled_madd_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, (buf_type)beta, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U << 1, num_bits_B);
                        // normal format (64bit source)
                        } else {
                            scaled_unpack<buf_type, base_type_B>
                                    (b_buf, B.packed_data, (buf_type)beta, B.get_loc(k, j), 
                                    k_inc, j_inc, M_U, num_bits_B);
                        }
                    }
                    
                    // execute mini-mmm on buffers with unpacked data
                    mini_mmm<buf_type, res_buf_type, mmm_version>(a_buf, b_buf, c_buf, 
                            i_inc, j_inc, k_inc, M_U, P_U);
                }
                
                // pack result to D
                pack<res_buf_type, base_type_D>
                        (c_buf, D.packed_data, D.get_loc(i, j), 
                        i_inc, j_inc, M_U, res_bits);
            }
        }
        
        // used buffers
        free(a_buf);
        free(b_buf);
        free(c_buf);

        // reset tiles after use
        if (use_amx) _tile_release();
    }

    return D;
}
