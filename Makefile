CC=g++
AS=as
SRC_DIR=src
ODIR=build
CFLAGS=-march=native -O3 -Iinclude
ASM_FILES=mini_mmm_asm_8_32.s mini_mmm_asm_16_32.s mini_mmm_asm_64_64.s mini_mmm_asm_64_64_ifma.s
CPP_FILES=APIB_packing.cpp APIB_unpacking.cpp APIB_scaled_unpacking.cpp APIB_mmm.cpp
ASM_PATHS=$(addprefix $(SRC_DIR)/,$(ASM_FILES))
CPP_PATHS=$(addprefix $(SRC_DIR)/,$(CPP_FILES))
APIB_OBJS=$(patsubst $(SRC_DIR)/%.cpp,$(ODIR)/%.o,$(CPP_PATHS)) $(patsubst $(SRC_DIR)/%.s,$(ODIR)/%.o,$(ASM_PATHS))

NUM_THREADS?=1
MAX_NUM_THREADS?=12
TASK_MASK?=0x1
MAX_BLOCK_SIZE?=256

DO_APIB?=1
DO_APIB_Slice?=0
DO_EIGEN?=0
DO_MKL?=0
DO_LUT?=0

ifeq ($(NUM_THREADS), 1)
	THREADS=SEQ
else
	THREADS=PAR
endif

BASE_FLAGS=-march=native -O3 -fno-strict-aliasing -fno-tree-vectorize -I/${PAPI_DIR}/include -L/${PAPI_DIR}/lib -lpapi -Iinclude

APIB_FLAGS_SEQ= $(BASE_FLAGS) -Llib -lAPIB
APIB_FLAGS_PAR= $(APIB_FLAGS_SEQ) -fopenmp
APIB_FLAGS=$(APIB_FLAGS_$(THREADS))

EIGEN_FLAGS_SEQ=$(BASE_FLAGS) -I/${EIGEN_DIR}
EIGEN_FLAGS_PAR=$(EIGEN_FLAGS_SEQ) -fopenmp
EIGEN_FLAGS=$(EIGEN_FLAGS_$(THREADS))

MKL_FLAGS_SEQ=-m64  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -I/${PAPI_DIR}/include -L/${PAPI_DIR}/lib -lpapi -Iinclude
MKL_FLAGS_PAR=-m64 -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -I/${PAPI_DIR}/include -L/${PAPI_DIR}/lib -lpapi -Iinclud -Iinclude
MKL_FLAGS=$(MKL_FLAGS_$(THREADS))

$(ODIR):
	mkdir -p $(ODIR)

$(ODIR)/%.o: $(SRC_DIR)/%.cpp | $(ODIR)
	$(CC) $(BASE_FLAGS) -c $< -o $@

$(ODIR)/%.o: $(SRC_DIR)/%.s | $(ODIR)
	$(AS) $< -o $@

lib:
	mkdir -p lib

lib/libAPIB.a: $(APIB_OBJS) | lib
	ar rcs $@ $^

init: lib/libAPIB.a
	mkdir -p logs
	mkdir -p results/rooflines
	mkdir -p results/measurements/lut
	mkdir -p results/measurements/square
	mkdir -p results/measurements/thin_1
	mkdir -p results/measurements/thin_2
	mkdir -p results/measurements/slices

reproduce_results_all:
	make DO_APIB=1 DO_APIB_Slice=1 DO_EIGEN=1 DO_MKL=1 DO_LUT=1 reproduce_results

reproduce_results_apib:
	make DO_APIB=1 DO_APIB_Slice=0 DO_EIGEN=0 DO_MKL=0 DO_LUT=0 reproduce_results

reproduce_results_apib_slice:
	make DO_APIB=0 DO_APIB_Slice=1 DO_EIGEN=0 DO_MKL=0 DO_LUT=0 reproduce_results

reproduce_results_eigen:
	make DO_APIB=0 DO_APIB_Slice=0 DO_EIGEN=1 DO_MKL=0 DO_LUT=0 reproduce_results

reproduce_results_mkl:
	make DO_APIB=0 DO_APIB_Slice=0 DO_EIGEN=0 DO_MKL=1 DO_LUT=0 reproduce_results

reproduce_results_lut:
	make NUM_THREADS=1 DO_APIB=0 DO_APIB_Slice=0 DO_EIGEN=0 DO_MKL=0 DO_LUT=1 bench_perf

reproduce_results:
	make NUM_THREADS=1 DO_APIB=$(DO_APIB) DO_APIB_Slice=$(DO_APIB_Slice) DO_EIGEN=$(DO_EIGEN) DO_MKL=$(DO_MKL) DO_LUT=$(DO_LUT) bench_perf
	make NUM_THREADS=4 DO_APIB=$(DO_APIB) DO_APIB_Slice=$(DO_APIB_Slice) DO_EIGEN=$(DO_EIGEN) DO_MKL=$(DO_MKL) DO_LUT=$(DO_LUT) bench_perf
	make NUM_THREADS=8 DO_APIB=$(DO_APIB) DO_APIB_Slice=$(DO_APIB_Slice) DO_EIGEN=$(DO_EIGEN) DO_MKL=$(DO_MKL) DO_LUT=$(DO_LUT) bench_perf
	make NUM_THREADS=12 DO_APIB=$(DO_APIB) DO_APIB_Slice=$(DO_APIB_Slice) DO_EIGEN=$(DO_EIGEN) DO_MKL=$(DO_MKL) DO_LUT=$(DO_LUT) bench_perf

bench_blocks:
	$(CC) src/benchmark.cpp -o src/benchmark $(CFLAGS)
	./src/benchmark 1 $(MAX_BLOCK_SIZE)
	rm -f src/res.txt src/benchmark

bench_bandwidth:
	$(CC) src/benchmark.cpp -o src/benchmark $(CFLAGS)
	./src/benchmark 2 $(MAX_NUM_THREADS)
	rm -f src/res.txt src/benchmark

bench_perf:
	$(CC) src/benchmark.cpp -o src/benchmark $(CFLAGS)
	./src/benchmark 0 $(NUM_THREADS) $(DO_APIB) $(DO_APIB_Slice) $(DO_EIGEN) $(DO_MKL) $(DO_LUT)
	rm -f src/res.txt src/benchmark

gen_apib:
	$(CC) src/gen.cpp -o src/gen $(APIB_FLAGS)
	OMP_NUM_THREADS=$(NUM_THREADS) taskset $(TASK_MASK) ./src/gen
	rm -f src/gen src/gen.cpp

gen_eigen:
	$(CC) src/gen.cpp -o src/gen $(EIGEN_FLAGS)
	OMP_NUM_THREADS=$(NUM_THREADS) taskset $(TASK_MASK) ./src/gen
	rm -f src/gen src/gen.cpp

gen_mkl:
	$(CC) src/gen.cpp -o src/gen $(MKL_FLAGS)
	OMP_NUM_THREADS=$(NUM_THREADS) taskset $(TASK_MASK) ./src/gen
	rm -f src/gen src/gen.cpp
	
test:
	$(CC) tests/test.cpp -o tests/test $(APIB_FLAGS_SEQ)
	./tests/test
	$(CC) tests/test.cpp -o tests/test $(APIB_FLAGS_PAR)
	OMP_NUM_THREADS=4 ./tests/test
	rm -f tests/test

bench_amx:
	$(CC) util/bench_amx.cpp -o util/bench_amx util/func_amx.s $(APIB_FLAGS_SEQ)
	./util/bench_amx
	rm util/bench_amx

clean:
	rm -f lib/* build/* logs/*
	
.PHONY: clean
