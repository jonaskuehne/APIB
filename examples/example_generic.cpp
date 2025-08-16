/**
 * @file example_generic.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Shows how to use the library APIB with generic matrices
 * 
 */

#include <cstdint>
#include <cstdlib>
#include "APIB/APIB.h"

int main() {
    // specifications for matrix
    const int n = 100;
    const int m = 200;
    const int p = 300;

    const int num_bits_A = 5;
    const int num_bits_B = 6;
    const int num_bits_C = 9;

    const int num_bits_alpha = 3;
    const int num_bits_beta = 2;

    // create buffers
    uint8_t* a_buf = (uint8_t*)malloc(n*p*sizeof(uint8_t));
    uint8_t* b_buf = (uint8_t*)malloc(p*m*sizeof(uint8_t));
    uint16_t* c_buf = (uint16_t*)malloc(n*m*sizeof(uint16_t));
    uint32_t* d_buf = (uint32_t*)malloc(n*m*sizeof(uint32_t));

    // initialize them in some meaningful way...
    
    // create matrices
    APIB_Mat<num_bits_A, n, p, uint8_t> mat_A(a_buf);
    APIB_Mat<num_bits_B, p, m, uint8_t> mat_B(b_buf);
    APIB_Mat<num_bits_C, n, m, uint16_t> mat_C(c_buf);

    // perform mmm
    uint8_t alpha = 4;
    uint8_t beta = 3;
    auto mat_D = mmm<num_bits_alpha, num_bits_beta>(mat_A, mat_B, mat_C, alpha, beta);

    // extract result
    mat_D.template unpack_data<uint32_t>(d_buf, m);

    // use d_buf for your computation ...
    
    // free memory
    free(a_buf);
    free(b_buf);
    free(c_buf);
    free(d_buf);
}
