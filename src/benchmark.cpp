/**
 * @file benchmark.cpp
 * @author Jonas Kühne (jonas.kuehne@proton.me)
 * @brief Main file to run all benchmarks from
 * 
 */

#include <algorithm>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "APIB/APIB_consts.h"

#define SLICE_DIM 32
#define APIB_SKELETON "src/skeletons/APIB_bench_skeleton.cpp"
#define APIB_SLICE_SKELETON "src/skeletons/APIB_slice_bench_skeleton.cpp"
#define EIGEN_SKELETON "src/skeletons/eigen_bench_skeleton.cpp"
#define MKL_SKELETON "src/skeletons/mkl_bench_skeleton.cpp"
#define BANDWIDTH_SKELETON "src/skeletons/bandwidth_measure_bench_skeleton.cpp"
#define LUT_SKELETON "src/skeletons/LUT_bench_skeleton.cpp"
#define CONST_FILE "include/APIB/APIB_consts.h"
#define LOG_DIR "logs/"
#define MEASURE_DIR "results/measurements/"

enum Bits {b8_32, b16_32, b64_64, b64_64_ifma, b8_32_lut};
enum Skel {APIB, APIB_Slice, Eigen, MKL, Bandwidth, LUT};
enum Experiment {square, thin_1, thin_2, slices};

// computes taskmask for number of threads
int get_task_mask(int num_threads) {
    return (1 << num_threads) - 1;
}

// reads performance and operational intensity from a res file
void read_perf(const char* res_name, double* store) {
    FILE* res = fopen(res_name, "r");
    int ret = (fscanf(res, "%lf\n", store) == EOF);
    ret |= (fscanf(res, "%lf\n", store + 1) == EOF);
    fclose(res);

    if (ret) {
        printf("read res_file went wrong\n");
        exit(1);
    }
}

// generates gen file from skeleton and replaces necessary texts
void gen(const char* res_name, const char* skeleton_name, const char* gen_name, bool verbose, Bits bits, int* params, int* bits_slices) {
    int buf_size = 200;
    char buf[buf_size];
    FILE* skeleton = fopen(skeleton_name, "r");
    FILE* gen = fopen(gen_name, "w");

    int include = 1;
    int num_slices = ceil((double)params[0] / SLICE_DIM);
    while(fgets(buf, buf_size, skeleton) != NULL) {
        std::string str (buf);
        if (str.find("/*GEN INSERT RESFILE*/") != std::string::npos) {
            fprintf(gen, "\"%s\"", res_name);
        } else if (str.find("/*GEN MAT_N*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[0]);
        } else if (str.find("/*GEN MAT_M*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[1]);
        } else if (str.find("/*GEN MAT_P*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[2]);
        } else if (str.find("/*GEN BITS_A*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[3]);
        } else if (str.find("/*GEN BITS_B*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[4]);
        } else if (str.find("/*GEN BITS_C*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[5]);
        } else if (str.find("/*GEN BITS_ALPHA*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[6]);
        } else if (str.find("/*GEN BITS_BETA*/") != std::string::npos) {
            fprintf(gen, "%d\n", params[7]);
        // type of eigen matrices
        } else if (str.find("/*GEN SCALAR_TYPE*/") != std::string::npos) {
            if (bits == b8_32 || bits == b16_32) {
                fprintf(gen, "uint32_t\n");
            } else {
                fprintf(gen, "uint64_t\n");
            }
        // type of A and B in mkl
        } else if (str.find("/*GEN OP_TYPE*/") != std::string::npos) {
            if (bits == b8_32) {
                fprintf(gen, "uint8_t\n");
            } else {
                fprintf(gen, "uint16_t\n");
            }
        // function to use in mkl
        } else if (str.find("/*GEN MKL_MMM*/") != std::string::npos) {
            if (bits == b8_32) {
                fprintf(gen, "cblas_gemm_s8u8s32(CblasRowMajor, CblasNoTrans, CblasNoTrans, CblasFixOffset, n, m, p, (float)alpha, A, p, 0, B, m, 0, (float)beta, (int *)C, m, &l);\n");
            } else {
                fprintf(gen, "cblas_gemm_s16s16s32(CblasRowMajor, CblasNoTrans, CblasNoTrans, CblasFixOffset, n, m, p, (float)alpha, (short *)A, p, 0, (short *)B, m, 0, (float)beta, (int *)C, m, &l);\n");
            }
        // print output
        } else if (!verbose && str.find("/*GEN VERBOSE START*/") != std::string::npos) {
            include = 0;
        } else if (!verbose && str.find("/*GEN VERBOSE END*/") != std::string::npos) {
            include = 1;
        } else if (str.find("/*GEN N_SLICES*/") != std::string::npos) {
            for (int i = 0; i < num_slices; i++) {
                fprintf(gen, "\tconst int n_%d = %d;\n", i, std::min(SLICE_DIM, params[0] - i*SLICE_DIM));
            }
        } else if (str.find("/*GEN BITS_SLICES*/") != std::string::npos) {
            if (bits_slices != NULL) {
                for (int i = 0; i < num_slices; i++) {
                    fprintf(gen, "\tconst int bits_A_%d = %d;\n", i, bits_slices[i]);
                } 
            } else {
                int max_bits = params[3];
                int min_bits = std::max(1, max_bits / 2);
                int curr_bits = max_bits;
                for (int i = 0; i < num_slices; i++) {
                    fprintf(gen, "\tconst int bits_A_%d = %d;\n", i, curr_bits);
                    if (curr_bits <= min_bits) {
                        curr_bits = max_bits;
                    } else {
                        curr_bits--;
                    }
                }
            }
        } else if (str.find("/*GEN CREATE_SPLIT_BUF*/") != std::string::npos) {
            for (int i = 0; i < num_slices; i++) {
                fprintf(gen, "uint64_t* A_data_%d = (uint64_t*)malloc(n_%d*p*sizeof(uint64_t));\n", i, i);
                fprintf(gen, "fill_matrix<uint64_t>(A_data_%d, n_%d, p, bits_A_%d);\n", i, i, i);
            }
        } else if (str.find("/*GEN CREATE_SLICES*/") != std::string::npos) {
            for (int i = 0; i < num_slices; i++) {
                fprintf(gen, "using S%d = APIB_Slice_Spec<bits_A_%d, n_%d, uint64_t>;\n", i, i, i);
            }
        } else if (str.find("/*GEN CREATE_SPLIT_MAT_OBJ*/") != std::string::npos) {
            fprintf(gen, "APIB_Slice_Mat<p");
            for (int i = 0; i < num_slices; i++) {
                fprintf(gen, ", S%d", i);
            }
            fprintf(gen, "> mat_A(");

            for (int i = 0; i < num_slices; i++) {
                if (i > 0) {
                    fprintf(gen, ", ");
                }
                fprintf(gen, "A_data_%d", i);
            }
            fprintf(gen, ");\n");

        } else if (str.find("/*GEN FREE_SPLIT_BUF*/") != std::string::npos) {
            for (int i = 0; i < num_slices; i++) {
                fprintf(gen, "free(A_data_%d);\n", i);
            }
        } else if (str.find("/*GEN CREATE_SPLIT_MAT_LIST*/") != std::string::npos) {
            fprintf(gen, "APIB_Slice_Mat<p");
            for (int i = 0; i < num_slices; i++) {
                fprintf(gen, ", S%d", i);
            }
            fprintf(gen, "> A[limit];\n");
        } else if (include && str.find("/*GEN") == std::string::npos) {
            fputs(buf, gen);
        }
    }

    fclose(skeleton);
    fclose(gen);
}

// executes a run of the config to be benchmarked
void run(const char* res_name, const char* skeleton_name, bool verbose, Bits bits, Skel skel, int* params, int taskmask, int num_threads, int* bits_slices) {
    gen(res_name, skeleton_name, "src/gen.cpp", verbose, bits, params, bits_slices);
    int ret;
    
    // construct command
    std::string s("make TASK_MASK=");
    s += std::to_string(taskmask) + " NUM_THREADS="+std::to_string(num_threads) + " ";

    switch (skel) {
        case APIB:
	    case APIB_Slice:
	    case LUT:
            ret = system((s + "gen_apib").c_str());
            break;
        case Eigen:
            ret = system((s + "gen_eigen").c_str());
            break;
        case Bandwidth:
        case MKL:
            ret = system((s + "gen_mkl").c_str());
            break;
        default:
            printf("skel not implemented\n");
            exit(1);
    }

    if (ret) {
        printf("make went wrong\n");
        exit(1);
    }

}

// example bit sizes per buffer type
int get_params(int* params, Bits bits) {
    int p_max;

    switch (bits) {
        case b8_32:
            params[3] = 2;
            params[4] = 2;
            params[5] = 2;
            params[6] = 2;
            params[7] = 2;
            p_max = 1 << 25;
            break;
        case b16_32:
            params[3] = 9;
            params[4] = 8;
            params[5] = 10;
            params[6] = 1;
            params[7] = 1;
            p_max = 1 << 13;
            break;
        case b64_64:
            params[3] = 18;
            params[4] = 33;
            params[5] = 20;
            params[6] = 2;
            params[7] = 2;
            p_max = 1 << 10;
            break;
        case b64_64_ifma:
            params[3] = 3;
            params[4] = 33;
            params[5] = 20;
            params[6] = 2;
            params[7] = 2;
            p_max = 1 << 13;
            break;
        case b8_32_lut:
            params[3] = 1;
            params[4] = 1;
            params[5] = 1;
            params[6] = 1;
            params[7] = 1;
            p_max = 1 << 25;
            break;
        default:
            printf("op/res bit combination not implemented\n");
            exit(1);

    }

    return p_max;
}

// writes constants to CONST_FILE
void write_params(int* block_8_32, int* block_16_32, int* block_64_64_ifma, int* block_64_64) {
    FILE* param_file = fopen(CONST_FILE, "w");
    // print header
    fprintf(param_file, "/**\n * @file APIB_consts.h\n * @author Jonas Kühne (jonas.kuehne@proton.me)\n * @brief Contains benchmarked parameters for multiplication of APIB matrix objects\n * \n */\n\n#pragma once\n\n"); 
    // get max of M_U, P_U for block width
    int block_max[] = {std::max(block_8_32[1], block_8_32[2]), 
        std::max(block_16_32[1], block_16_32[2]), 
        std::max(block_64_64_ifma[1], block_64_64_ifma[2]), 
        std::max(block_64_64[1], block_64_64[2])};

    int block_width = std::max(std::max(block_max[0], block_max[1]), std::max(block_max[2], block_max[3]));

    int block_exp = ((int)ceil(log2(block_width)));

    fprintf(param_file, "#define APIB_BLOCK_EXP %d\n\n", block_exp);

    fprintf(param_file, "#define APIB_N_U_8_32 %d\n#define APIB_M_U_8_32 %d\n#define APIB_P_U_8_32 %d\n\n", block_8_32[0], block_8_32[1], block_8_32[2]);
    fprintf(param_file, "#define APIB_N_U_16_32 %d\n#define APIB_M_U_16_32 %d\n#define APIB_P_U_16_32 %d\n\n", block_16_32[0], block_16_32[1], block_16_32[2]);
    fprintf(param_file, "#define APIB_N_U_64_64_IFMA %d\n#define APIB_M_U_64_64_IFMA %d\n#define APIB_P_U_64_64_IFMA %d\n\n", block_64_64_ifma[0], block_64_64_ifma[1], block_64_64_ifma[2]);
    fprintf(param_file, "#define APIB_N_U_64_64 %d\n#define APIB_M_U_64_64 %d\n#define APIB_P_U_64_64 %d\n\n", block_64_64[0], block_64_64[1], block_64_64[2]);
    fclose(param_file);
}

// runner to benchmark best block sizes
void run_bench_block_sizes(const char* store_name, Bits bits, int lim, int* dest) {
    int params[8];
    int p_max = get_params(params, bits);
    // size to be used
    params[0] = 2048;
    params[1] = 2048;
    params[2] = std::min(p_max, 2048);

    int N_U_fast = 0;
    int M_U_fast = 0;
    int P_U_fast = 0;
    double perf = 0;
    double op_int = 0;

    for (int N_U = 64; N_U <= 2*lim; N_U += 16) {
        for (int M_U = 128; M_U <= lim; M_U *= 2) {
            for (int P_U = 128; P_U <= lim; P_U *= 2) {

                int blocks[] = {N_U, M_U, P_U};
                write_params(blocks, blocks, blocks, blocks);
                run("src/res.txt", APIB_SKELETON, 0, bits, APIB, params, 0x1, 1, NULL); 
                double vals[2];
                read_perf("src/res.txt", vals);
                double p = vals[0];
                double intens = vals[1];

                if (p > perf) {
                    FILE* store_file = fopen(store_name, "a");
                    perf = p;
                    op_int = intens;
                    N_U_fast = N_U;
                    M_U_fast = M_U;
                    P_U_fast = P_U;
                    fprintf(store_file, "N_U=%d, M_U=:%d, P_U=%d, perf=%lf ops/cycle, op_int=%lf ops/byte\n", N_U, M_U, P_U, p, intens); 
                    fclose(store_file);
                }
            }
        }
    }

    dest[0] = N_U_fast;
    dest[1] = M_U_fast;
    dest[2] = P_U_fast;
}

// finds best block sizes N_U, M_U, P_U for all types by exhaustive search
void bench_block_sizes(int max_block_size) {
    int blocks_8_32[3];
    int blocks_16_32[3];
    int blocks_64_64_ifma[3];
    int blocks_64_64[3];

    std::string s(LOG_DIR);
    run_bench_block_sizes((s + "8_32_bench_param.txt").c_str(), b8_32, max_block_size, blocks_8_32);
    run_bench_block_sizes((s + "16_32_bench_param.txt").c_str(), b16_32, max_block_size, blocks_16_32);
    run_bench_block_sizes((s + "64_64_IFMA_bench_param.txt").c_str(), b64_64_ifma, max_block_size, blocks_64_64_ifma);
    run_bench_block_sizes((s + "64_64_bench_param.txt").c_str(), b64_64, max_block_size, blocks_64_64);

    write_params(blocks_8_32, blocks_16_32, blocks_64_64_ifma, blocks_64_64);
}

// runner for experiments where dimensions are fix and bits in slices change
void run_slices(const char* store_name, Bits bits, Skel skel, int n, int m, int p, int taskmask, int num_threads) {
    char skel_name[100];
    switch (skel) {
        case APIB_Slice:
            strcpy(skel_name, APIB_SLICE_SKELETON);
            break;
        case MKL:
            strcpy(skel_name, MKL_SKELETON);
            break;
        case Eigen:
            strcpy(skel_name, EIGEN_SKELETON);
            break;
        default:
            printf("skel not implemented\n");
            exit(1);
    }

    FILE* store_file = fopen(store_name, "w");
    char res_name[] = "src/res.txt";
    int params[8];
    int p_max = get_params(params, bits);
    params[0] = n;
    params[1] = m;
    params[2] = p;

    if (skel == APIB_Slice) {

        int bits_A_1[] = {2, 2, 2, 2, 4, 6, 8};
        int bits_A_2[] = {2, 4, 6, 8, 8, 8, 8};

        int its = sizeof(bits_A_1) / sizeof(int);

        for (int i = 0; i < its; i++) {
            int bits_A[] = {bits_A_1[i], bits_A_2[i]}; 
            run(res_name, skel_name, 1, bits, skel, params, taskmask, num_threads, bits_A); 
            double vals[2];
            read_perf(res_name, vals);
            double perf = vals[0];
            double bandwidth = vals[1];
            fprintf(store_file, "(%d, %d), %d, %d, %d\n", bits_A[0], bits_A[1], params[0], params[1], params[2]);
            fprintf(store_file, "%lf\n", perf);
            fprintf(store_file, "%lf\n", bandwidth);
        }
        fclose(store_file);
    } else {
        run(res_name, skel_name, 1, bits, skel, params, taskmask, num_threads, NULL); 
        double vals[2];
        read_perf(res_name, vals);
        double perf = vals[0];
        double bandwidth = vals[1];
        fprintf(store_file, "%d, %d, %d\n", params[0], params[1], params[2]);
        fprintf(store_file, "%lf\n", perf);
        fprintf(store_file, "%lf\n", bandwidth);
        fclose(store_file);
    }
}


// runner to benchmark the different implementation on different matrix sizes and shapes
void run_bench_mat_sizes(const char* store_name, Bits bits, Skel skel, Experiment exper, bool ignore_p_max, int min, int max, int inc, int taskmask, int num_threads) {
    char skel_name[100];
    switch (skel) {
        case APIB:
            strcpy(skel_name, APIB_SKELETON);
            break;
        case APIB_Slice:
            strcpy(skel_name, APIB_SLICE_SKELETON);
            break;
        case Eigen:
            strcpy(skel_name, EIGEN_SKELETON);
            break;
        case MKL:
            strcpy(skel_name, MKL_SKELETON);
            break;
        case Bandwidth:
            strcpy(skel_name, BANDWIDTH_SKELETON);
            break;
        case LUT:
            strcpy(skel_name, LUT_SKELETON);
            break;
        default:
            printf("skel not implemented\n");
            exit(1);
    }

    int params[8];
    int p_max = get_params(params, bits);
    
    FILE* store_file = fopen(store_name, "w");

    char res_name[] = "src/res.txt";
    for (int s = min; s <= max;/* see end of body of loop */) {
        switch (exper) {
            case square:
                params[0] = s;
                params[1] = s;
                if (ignore_p_max) {
                    params[2] = s;
                } else {
                    params[2] = std::min(s, p_max);
                }
                break;
            case thin_1:
                params[0] = s;
                params[1] = max;
                if (ignore_p_max) {
                    params[2] = max;
                } else {
                    params[2] = std::min(max, p_max);
                }
                break;
            case thin_2:
                params[0] = 8;
                params[1] = s;
                if (ignore_p_max) {
                    params[2] = s;
                } else {
                    params[2] = std::min(s, p_max);
                }
                break;
            default:
                printf("experiment not implemented\n");
                exit(1);
        }

        run(res_name, skel_name, 1, bits, skel, params, taskmask, num_threads, NULL); 
        double vals[2];
        read_perf(res_name, vals);
        double perf = vals[0];
        double bandwidth = vals[1];
        fprintf(store_file, "%d, %d, %d\n", params[0], params[1], params[2]);
        fprintf(store_file, "%lf\n", perf);
        fprintf(store_file, "%lf\n", bandwidth);
        printf("Done with s=%d\n", s);

        switch (exper) {
            case square:
            case thin_2:
                s+=inc;
                break;
            case thin_1:
                s*=inc;
                break;
            default:
                printf("experiment not implemented\n");
                exit(1);
        }
    }

    fclose(store_file);
}

// generates arguments function run_bench_mat_sizes and executes it
void exec_bench_mat(Skel skel, Bits bit, int threads, bool for_lut) {
    std::string s(MEASURE_DIR);
    std::string square_s = s + "square/";
    std::string thin_s_1 = s + "thin_1/";
    std::string thin_s_2 = s + "thin_2/";
    std::string slices_s = s + "slices/";
    std::string lut_s = s + "lut/";

    // convert skeleton to string
    std::string skel_str;
    switch (skel) {
        case APIB:
            skel_str = std::string("APIB");
            break;
        case APIB_Slice:
            skel_str = std::string("APIB_Sliced");
            break;
        case Eigen:
            skel_str = std::string("Eigen");
            break;
        case MKL:
            skel_str = std::string("MKL");
            break;
        case LUT:
            skel_str = std::string("LUT");
            break;
        default:
            // could be Bandwidth, should be ignored here
            break;
    }

    // convert bit combination to string
    std::string bit_str;
    bool ignore_max = 0;
    switch (bit) {
        case b8_32:
        case b8_32_lut:
            bit_str = std::string("8_32");
            break;
        case b16_32:
            bit_str = std::string("16_32");
            break;
        case b64_64_ifma:
            bit_str = std::string("64_64_IFMA");
            break;
        case b64_64:
            bit_str = std::string("64_64");
            ignore_max = 1;
            break;
        default:
            printf("bit not implemented\n");
            exit(1);
    }

    // dimensions and bitmask to taskset
    int square_start;
    int square_inc;
    int square_max;
    int task_mask;
    int block_width = 1 << APIB_BLOCK_EXP;

    if (threads <= 1) {
        if (for_lut) {
            square_start = 1024;
        } else {
            square_start = 1000;
        }
        square_inc = square_start;
        square_max = 5 * square_start;
    } else {
        square_start = threads * block_width;
        square_inc = square_start;
        
        if (threads >= 8) {
            square_max = 4 * square_start;
        } else {
            square_max = 5 * square_start;
        }
    }
    task_mask = get_task_mask(threads);

    std::string name = skel_str + "_" + bit_str + "_" + std::to_string(threads) + "_threads.txt";
    if (for_lut) {
        run_bench_mat_sizes((lut_s + name).c_str(), bit, skel, square, ignore_max, square_start, square_max, square_inc, task_mask, threads);
    } else {
        run_bench_mat_sizes((square_s + name).c_str(), bit, skel, square, ignore_max, square_start, square_max, square_inc, task_mask, threads);
        run_bench_mat_sizes((thin_s_1 + name).c_str(), bit, skel, thin_1, ignore_max, 3, 3072, 2, task_mask, threads);
        run_bench_mat_sizes((thin_s_2 + name).c_str(), bit, skel, thin_2, ignore_max, 1024, 3072, 1024, task_mask, threads);
        if (bit == b8_32 && (skel == MKL || skel == APIB_Slice || skel == Eigen)) {
            run_slices((slices_s + name).c_str(), bit, skel, 64, 2048, 1, task_mask, threads);
        }
    }
}

// reproduces result for given number of threads and library selection
void bench_mat_sizes(int threads, bool do_apib, bool do_apib_slice, bool do_eigen, bool do_mkl, bool do_lut) {
    // all bit combinations
    Bits all_bits[] = {b8_32, b16_32, b64_64_ifma, b64_64};

    // run all for apib
    for (int bit_ind = 0; bit_ind < sizeof(all_bits) / sizeof(Bits); bit_ind++) {
        if (do_apib) exec_bench_mat(APIB, all_bits[bit_ind], threads, false);
        if (do_apib_slice) exec_bench_mat(APIB_Slice, all_bits[bit_ind], threads, false);
    } 

    // run those supported by mkl
    if (do_mkl) {
        exec_bench_mat(MKL, b8_32, threads, false);
        exec_bench_mat(MKL, b16_32, threads, false);
    }

    // run 8->32 and 64->64 ifma, 16->32 and 64->64 same, copies file at end
    if (do_eigen) {
        exec_bench_mat(Eigen, b8_32, threads, false);
        exec_bench_mat(Eigen, b64_64_ifma, threads, false);

        std::string s(MEASURE_DIR);
        std::string name_8_32 = std::string("Eigen_8_32_") + std::to_string(threads) + "_threads.txt";
        std::string name_64_64_ifma = std::string("Eigen_64_64_IFMA_") + std::to_string(threads) + "_threads.txt";

        std::string name_16_32 = std::string("Eigen_16_32_") + std::to_string(threads) + "_threads.txt";
        std::string name_64_64 = std::string("Eigen_64_64_") + std::to_string(threads) + "_threads.txt";

        std::string start("cp ");

        int ret = system((start + s + "square/" + name_8_32 + " " + s + "square/" + name_16_32).c_str());
        ret |= system((start + s + "thin_1/" + name_8_32 + " " + s + "thin_1/" + name_16_32).c_str());
        ret |= system((start + s + "thin_2/" + name_8_32 + " " + s + "thin_2/" + name_16_32).c_str());
        ret |= system((start + s + "square/" + name_64_64_ifma + " " + s + "square/" + name_64_64).c_str());
        ret |= system((start + s + "thin_1/" + name_64_64_ifma + " " + s + "thin_1/" + name_64_64).c_str());
        ret |= system((start + s + "thin_2/" + name_64_64_ifma + " " + s + "thin_2/" + name_64_64).c_str());

        if (ret) {
            printf("make went wrong\n");
            exit(1);
        }
    }

    // run LUT stuff, only for one thread
    if (do_lut && threads == 1) {
        exec_bench_mat(LUT, b8_32_lut, 1, true);
        exec_bench_mat(APIB, b8_32_lut, 1, true);
    }
}

// benchmarks the memory bandwidth of the current cpu
void bench_bandwidth(int max_threads) {
    int params[8];
    
    std::string store_name = std::string(MEASURE_DIR) + "bandwidth.txt";
    FILE* store_file = fopen(store_name.c_str(), "w");
    int n_inc = 2000;
    char res_name[] = "src/res.txt";

    for (int num_threads = 1; num_threads <= max_threads; num_threads++) {
        double bandwidth = 0;
        int taskmask = get_task_mask(num_threads);
        for (int n = n_inc; n <= 5*n_inc; n+=n_inc) {
		printf("%d threads, n=%d\n", num_threads, n);
            params[0] = n;
            run(res_name, BANDWIDTH_SKELETON, 1, b8_32, Bandwidth, params, taskmask, num_threads, NULL); 
            double vals[2];
            read_perf(res_name, vals);
            bandwidth = std::max(vals[0], bandwidth);
        }
        bandwidth *= num_threads;
        fprintf(store_file, "%lf\n", bandwidth);
    }

    fclose(store_file);
}

int main(int argc, char *argv[]) {
    int util = std::atoi(argv[1]);
    int arg = std::atoi(argv[2]);
    if (util == 0) {
        bench_mat_sizes(arg, std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]), std::atoi(argv[6]), std::atoi(argv[7])); 
    } else if (util == 1) {
        bench_block_sizes(arg);
    } else {
        bench_bandwidth(arg);
    }
    
    return 0;
}

