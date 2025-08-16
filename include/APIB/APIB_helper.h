/**
 * @file APIB_helper.h
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Contains multiple helper functions used in APIB project
 * 
 */

#pragma once

#include <cstdint>
#include <variant>
#include <vector>
#include <stdio.h>
#include <immintrin.h>
#include <type_traits>
#include <unistd.h>
#include <sys/syscall.h>
#include "APIB/APIB_consts.h"

// constants used for AMX
#define MAX_ROWS 16
#define MAX_COLSB 64
#define ARCH_GET_XCOMP_PERM     0x1022
#define ARCH_REQ_XCOMP_PERM     0x1023
#define XFEATURE_XTILECFG       17
#define XFEATURE_XTILEDATA      18

// macros for compile-time computation
#define MAX(A, B) ((A > B) ? (A) : (B))
#define MIN(A, B) ((A > B) ? (B) : (A))
// turn to 1 for debug info
#define CHECK_TEMPLATE 0
#define CHECK_MATRIX 0
#define CHECK_RES_BITS 0
#define CHECK_TYPES 0

// crashes program if types are too large
#define CRASH_IF_ERROR 0

/**
 * @brief Struct for AMX config
 * 
 */
typedef struct __tile_config {
    uint8_t palette_id;
    uint8_t start_row;
    uint8_t reserved_0[14];
    uint16_t colsb[16]; 
    uint8_t rows[16]; 
} __tilecfg;

/**
 * @brief Struct containing location of element
 * 
 */
struct Packed_Loc {
    int packed_num;
    int eff_cols;
    int bit_offset;

    /**
     * @brief Construct a new Packed_Loc object
     * 
     * @param packed_num Position in memory
     * @param eff_cols Effective column number
     * @param bit_offset Offset in packed format
     */
    Packed_Loc(int packed_num, int eff_cols, int bit_offset) : 
            packed_num(packed_num), eff_cols(eff_cols), bit_offset(bit_offset){}
};

/**
 * @brief Determines the type of the operand (where A and B get unpacked to) buffers
 * 
 * @tparam bits_op Number of bits the operands (with scaling) take
 * @tparam bits_res Number of bits of the result
 */
template<int bits_op, int bits_res>
struct get_op_buf_type {
    using type = typename std::conditional<(bits_res > 32 || bits_op > 16), uint64_t, 
          typename std::conditional<bits_op <= 8, uint8_t, 
          typename std::conditional<bits_op < 16, uint16_t, 
          uint64_t>::type>::type>::type;
};

/**
 * @brief Determines the type of the result (where C and gets unpacked to and D gets packed from) buffer
 * 
 * @tparam bits_op Number of bits the operands (with scaling) take
 * @tparam bits_res Number of bits of the result
 */
template<int bits_op, int bits_res>
struct get_res_buf_type {
    using type = typename std::conditional<(bits_res > 32 || bits_op >= 16), uint64_t, uint32_t>::type;
};

/**
 * @brief Determines the base type of the APIB matrix object
 * 
 * @tparam num_bits Number of bits the elements take
 */
template<int num_bits>
struct get_base_type {
    using type = typename std::conditional<num_bits <= 32, uint32_t, 
          uint64_t>::type;
};

/**
 * @brief Computes ceiled log2 at compile time
 * 
 * @tparam x Argument to log
 */
template<int x>
struct Ceil_Log2 {
    static const constexpr unsigned int value = (x <= 1) ? 0 : 1 + Ceil_Log2<(x + 1) / 2>::value;
};

/**
 * @brief Computes the number of bits the unpacked elements of A and B (including scaling) take
 * 
 * @tparam num_bits_beta Number of bits of b
 * @tparam num_bits_A Number of bits of elements of A
 * @tparam num_bits_B Number of bits of elements of B
 */
template<int num_bits_beta, int num_bits_A, int num_bits_B>
struct Min_Bits_Op {
    static const constexpr unsigned int value = MAX(MAX(num_bits_A, num_bits_B), num_bits_beta + MIN(num_bits_A, num_bits_B));
};

/**
 * @brief Computes the number of bits the elements of D <- a*C + b*AxB take
 * 
 * @tparam num_bits_alpha Number of bits of a
 * @tparam num_bits_beta Number of bits of b
 * @tparam num_bits_A Number of bits of elements of A
 * @tparam num_bits_B Number of bits of elements of B
 * @tparam num_bits_C Number of bits of elements of C
 * @tparam p Number of columns of A (or rows of B)
 */
template<int num_bits_alpha, int num_bits_beta, int num_bits_A, int num_bits_B, int num_bits_C, int p>
struct Min_Bits_Res {
    static const constexpr unsigned int value = 1 + MAX(num_bits_beta + num_bits_A + num_bits_B + Ceil_Log2<p>::value, num_bits_alpha + num_bits_C);
};

// checks if T is in Ts
template <typename T, typename... Ts>
struct is_contained : std::disjunction<std::is_same<T, Ts>...> {};

// appends unique type
// base
template <typename T, typename Tuple>
struct append_unique;
// recursion
template <typename T, typename... Ts>
struct append_unique<T, std::tuple<Ts...>> {
    using type = std::conditional_t<
        is_contained<T, Ts...>::value,
        std::tuple<Ts...>,
        std::tuple<Ts..., T>>;
};

/**
 * @brief Creates a tuple of the types Ts without duplicates
 * 
 * @tparam Ts List of types
 */
template <typename... Ts>
struct unique_types {
    using type = std::tuple<>;
};
// recursion
template <typename T, typename... Ts>
struct unique_types<T, Ts...> {
    using tail_unique = typename unique_types<Ts...>::type;
    using type = typename append_unique<T, tail_unique>::type;
};

/**
 * @brief Converts a std::tuple to a std::variant
 * 
 * @tparam Tuple Tuple to be converted
 */
template <typename Tuple>
struct tuple_to_variant;
// recursion
template <typename... Ts>
struct tuple_to_variant<std::tuple<Ts...>> {
    using type = std::variant<Ts...>;
};

/**
 * @brief Represents a slice of a matrix as a triple
 * 
 * @tparam num_bits Number of bits the elements in the slice take
 * @tparam num_rows Number of rows of the slice
 * @tparam slice_type Type of the buffer with which the slice will be initialized
 */
template <int num_bits, int num_rows, typename slice_type>
struct APIB_Slice_Spec {
    using Type = slice_type;
    static constexpr int bits = num_bits;
    static constexpr int rows = num_rows;
};

/**
 * @brief Computes the maximum bit of a list of APIB_Slice_Specs
 * 
 * @tparam specs List of APIB_Slice_Specs
 */
template <typename... specs>
struct get_max_bits;
// base
template <typename spec>
struct get_max_bits<spec> {
    static constexpr int value = spec::bits;
};
// recursion
template <typename First, typename Second, typename... Rest>
struct get_max_bits<First, Second, Rest...> {
    static constexpr int value = (First::bits > get_max_bits<Second, Rest...>::value)
                                 ? First::bits
                                 : get_max_bits<Second, Rest...>::value;
};

/**
 * @brief Computes the sum of rows of a list of APIB_Slice_Specs
 * 
 * @tparam Specs List of APIB_Slice_Specs
 */
template <typename... Specs>
struct get_total_rows;
// base
template <typename Spec>
struct get_total_rows<Spec> {
    static constexpr int value = Spec::rows;
};
// recursion
template <typename First, typename Second, typename... Rest>
struct get_total_rows<First, Second, Rest...> {
    static constexpr int value = First::rows + get_total_rows<Second, Rest...>::value;
};

/**
 * @brief Selects blocksize N_U from benchmarked parameters
 * 
 * @tparam num_bits_buf number of bits of operand buffer type
 * @tparam version 0 if can use IFMA for 64bit sources, 1 if cannot
 */
template<int num_bits_buf, int version>
struct Get_N_U {
    static const constexpr unsigned int value = 
        // 8bit, 32bit
        (num_bits_buf == 8) ? APIB_N_U_8_32 : 
        // 16bit, 32bit
        ((num_bits_buf == 16) ? APIB_N_U_16_32 : 
         // 64bit, 64bit
         (version ? APIB_N_U_64_64 : 
          // 64bit, 64bit ifma
          APIB_N_U_64_64_IFMA));
};

/**
 * @brief Selects blocksize M_U from benchmarked parameters
 * 
 * @tparam num_bits_buf number of bits of operand buffer type
 * @tparam version 0 if can use IFMA for 64bit sources, 1 if cannot
 */
template<int num_bits_buf, int version>
struct Get_M_U {
    static const constexpr unsigned int value = 
        // 8bit, 32bit
        (num_bits_buf == 8) ? APIB_M_U_8_32 : 
        // 16bit, 32bit
        ((num_bits_buf == 16) ? APIB_M_U_16_32 : 
         // 64bit, 64bit
         (version ? APIB_M_U_64_64 : 
          // 64bit, 64bit ifma
          APIB_M_U_64_64_IFMA));
};

/**
 * @brief Selects blocksize P_U from benchmarked parameters
 * 
 * @tparam num_bits_buf number of bits of operand buffer type
 * @tparam version 0 if can use IFMA for 64bit sources, 1 if cannot
 */
template<int num_bits_buf, int version>
struct Get_P_U {
    static const constexpr unsigned int value = 
        // 8bit, 32bit
        (num_bits_buf == 8) ? APIB_P_U_8_32 : 
        // 16bit, 32bit
        ((num_bits_buf == 16) ? APIB_P_U_16_32 : 
         // 64bit, 64bit
         (version ? APIB_P_U_64_64 : 
          // 64bit, 64bit ifma
          APIB_P_U_64_64_IFMA));
};

/**
 * @brief Initializes the tile config for AMX
 * 
 * @param tileinfo Config
 * @param rows Number of rows in tile
 * @param colsb Size of one row in bytes
 */
static void init_tile_config (__tilecfg *tileinfo, int rows, int colsb) {
    int i;
    tileinfo->palette_id = 1;
    tileinfo->start_row = 0;

    for (i = 0; i < 8; i++) {
        tileinfo->rows[i] =  rows;
        tileinfo->colsb[i] = colsb;
    }

    _tile_loadconfig (tileinfo);
}

/**
 * @brief Set the tiledata use object
 * 
 * @return true Worked
 * @return false Failed
 */
static bool set_tiledata_use() {
    // invoke syscall to set ARCH_SET_STATE_USE
    if (syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA)) {
        printf("\n Fail to do XFEATURE_XTILEDATA \n\n");
        return false;
    }
   return true;
}
