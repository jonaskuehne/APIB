# file mini_mmm_asm_64_64.s
# author Jonas Kühne (jonas.kuehne@proton.me)
# Contains kernels for 64->64

# args:         rdi, rsi, rdx, rcx, r8, r9
# return:       rax
# caller saved: rax, rdi, rsi, rdx, rcx, r8, r9, r10, and r11;
# callee saved: rbx, rsp, rbp, r12, r13, r14, and r15;
.text
.globl mini_mmm_asm_64_64

.globl mini_mmm_asm_64_64_clean_4
.globl mini_mmm_asm_64_64_clean_4_unr_16
.globl mini_mmm_asm_64_64_clean_4_unr_32

.globl mini_mmm_asm_64_64_clean_2
.globl mini_mmm_asm_64_64_clean_2_unr_16
.globl mini_mmm_asm_64_64_clean_2_unr_32

.globl mini_mmm_asm_64_64_clean_1
.globl mini_mmm_asm_64_64_clean_1_unr_16
.globl mini_mmm_asm_64_64_clean_1_unr_32
.globl mini_mmm_asm_64_64_clean_1_unr_64

# %rdi: a_buf + i1*P_U
# %rsi: b_buf + j1
# %rdx: c_buf + i1*M_U + j1
# %rcx: k_inc
# %r8: P_U
# %r9: M_U

# scratch: %r10, %r11, %rax

# 8x8 block
mini_mmm_asm_64_64:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 8), %zmm1
	vmovdqu32 (%r10, %r9, 8), %zmm2
    
    leaq (%rax, %r9, 8), %r10

	vmovdqu32 (%r11, %r9, 8), %zmm3

	leaq (%r10, %r9, 8), %r11

	vmovdqu32 (%rax, %r9, 8), %zmm4

	leaq (%r11, %r9, 8), %rax

	vmovdqu32 (%r10, %r9, 8), %zmm5
	vmovdqu32 (%r11, %r9, 8), %zmm6
	vmovdqu32 (%rax, %r9, 8), %zmm7

	jmp .loop_cond
	
	.loop_body:
		# load b vecs
		vmovdqu32 (%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 8), %r10
		leaq (%r10, %r8, 8), %r11
		leaq (%r11, %r8, 8), %rax

	    vpbroadcastq (%rdi), %zmm8
	    vpbroadcastq (%rdi, %r8, 8), %zmm9
	    vpbroadcastq (%r10, %r8, 8), %zmm10
        
        leaq (%rax, %r8, 8), %r10

	    vpbroadcastq (%r11, %r8, 8), %zmm11

	    leaq (%r10, %r8, 8), %r11

	    vpbroadcastq (%rax, %r8, 8), %zmm12

	    leaq (%r11, %r8, 8), %rax

	    vpbroadcastq (%r10, %r8, 8), %zmm13
	    vpbroadcastq (%r11, %r8, 8), %zmm14
	    vpbroadcastq (%rax, %r8, 8), %zmm15

        # mul
        vpmullq %zmm8, %zmm31, %zmm16
        vpmullq %zmm9, %zmm31, %zmm17
        vpmullq %zmm10, %zmm31, %zmm18
        vpmullq %zmm11, %zmm31, %zmm19
        vpmullq %zmm12, %zmm31, %zmm20
        vpmullq %zmm13, %zmm31, %zmm21
        vpmullq %zmm14, %zmm31, %zmm22
        vpmullq %zmm15, %zmm31, %zmm23

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1
		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3
		vpaddq %zmm4, %zmm20, %zmm4
		vpaddq %zmm5, %zmm21, %zmm5
		vpaddq %zmm6, %zmm22, %zmm6
		vpaddq %zmm7, %zmm23, %zmm7

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body

	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 8)
	vmovdqu32 %zmm2, (%r10, %r9, 8)
    
    leaq (%rax, %r9, 8), %r10

	vmovdqu32 %zmm3, (%r11, %r9, 8)

	leaq (%r10, %r9, 8), %r11

	vmovdqu32 %zmm4, (%rax, %r9, 8)

	leaq (%r11, %r9, 8), %rax

	vmovdqu32 %zmm5, (%r10, %r9, 8)
	vmovdqu32 %zmm6, (%r11, %r9, 8)
	vmovdqu32 %zmm7, (%rax, %r9, 8)

	ret

# 4x8 block
mini_mmm_asm_64_64_clean_4:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 8), %zmm1
	vmovdqu32 (%r10, %r9, 8), %zmm2
	vmovdqu32 (%r11, %r9, 8), %zmm3

	jmp .loop_cond_clean_4
	
	.loop_body_clean_4:
		# load b vecs
		vmovdqu32 (%rsi), %zmm30
		
		# load a vecs
		leaq (%rdi, %r8, 8), %r10
		leaq (%r10, %r8, 8), %r11

	    vpbroadcastq (%rdi), %zmm8
	    vpbroadcastq (%rdi, %r8, 8), %zmm9
	    vpbroadcastq (%r10, %r8, 8), %zmm10
	    vpbroadcastq (%r11, %r8, 8), %zmm11

        # mul
        vpmullq %zmm8, %zmm30, %zmm16
        vpmullq %zmm9, %zmm30, %zmm17
        vpmullq %zmm10, %zmm30, %zmm18
        vpmullq %zmm11, %zmm30, %zmm19

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1
		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_4:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_4

	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 8)
	vmovdqu32 %zmm2, (%r10, %r9, 8)
	vmovdqu32 %zmm3, (%r11, %r9, 8)

	ret

# 4x16 block
mini_mmm_asm_64_64_clean_4_unr_16:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 8), %zmm1
	vmovdqu32 (%r10, %r9, 8), %zmm2
	vmovdqu32 (%r11, %r9, 8), %zmm3

	vmovdqu32 64(%rdx), %zmm4
	vmovdqu32 64(%rdx, %r9, 8), %zmm5
	vmovdqu32 64(%r10, %r9, 8), %zmm6
	vmovdqu32 64(%r11, %r9, 8), %zmm7

	jmp .loop_cond_clean_4_unr_16
	
	.loop_body_clean_4_unr_16:
		# load b vecs
		vmovdqu32 (%rsi), %zmm30
		vmovdqu32 64(%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 8), %r10
		leaq (%r10, %r8, 8), %r11

	    vpbroadcastq (%rdi), %zmm8
	    vpbroadcastq (%rdi, %r8, 8), %zmm9
	    vpbroadcastq (%r10, %r8, 8), %zmm10
	    vpbroadcastq (%r11, %r8, 8), %zmm11

        # mul
        vpmullq %zmm8, %zmm30, %zmm16
        vpmullq %zmm9, %zmm30, %zmm17
        vpmullq %zmm10, %zmm30, %zmm18
        vpmullq %zmm11, %zmm30, %zmm19

        vpmullq %zmm8, %zmm31, %zmm20
        vpmullq %zmm9, %zmm31, %zmm21
        vpmullq %zmm10, %zmm31, %zmm22
        vpmullq %zmm11, %zmm31, %zmm23

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1
		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3

		vpaddq %zmm4, %zmm20, %zmm4
		vpaddq %zmm5, %zmm21, %zmm5
		vpaddq %zmm6, %zmm22, %zmm6
		vpaddq %zmm7, %zmm23, %zmm7

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_4_unr_16:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_4_unr_16

	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 8)
	vmovdqu32 %zmm2, (%r10, %r9, 8)
	vmovdqu32 %zmm3, (%r11, %r9, 8)

	vmovdqu32 %zmm4, 64(%rdx)
	vmovdqu32 %zmm5, 64(%rdx, %r9, 8)
	vmovdqu32 %zmm6, 64(%r10, %r9, 8)
	vmovdqu32 %zmm7, 64(%r11, %r9, 8)

	ret

# 4x32 block
mini_mmm_asm_64_64_clean_4_unr_32:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 8), %zmm1
	vmovdqu32 (%r10, %r9, 8), %zmm2
	vmovdqu32 (%r11, %r9, 8), %zmm3

	vmovdqu32 64(%rdx), %zmm4
	vmovdqu32 64(%rdx, %r9, 8), %zmm5
	vmovdqu32 64(%r10, %r9, 8), %zmm6
	vmovdqu32 64(%r11, %r9, 8), %zmm7

	vmovdqu32 128(%rdx), %zmm8
	vmovdqu32 128(%rdx, %r9, 8), %zmm9
	vmovdqu32 128(%r10, %r9, 8), %zmm10
	vmovdqu32 128(%r11, %r9, 8), %zmm11

	vmovdqu32 192(%rdx), %zmm12
	vmovdqu32 192(%rdx, %r9, 8), %zmm13
	vmovdqu32 192(%r10, %r9, 8), %zmm14
	vmovdqu32 192(%r11, %r9, 8), %zmm15

	jmp .loop_cond_clean_4_unr_32
	
	.loop_body_clean_4_unr_32:
		# load b vecs
		vmovdqu32 (%rsi), %zmm28
		vmovdqu32 64(%rsi), %zmm29
		vmovdqu32 128(%rsi), %zmm30
		vmovdqu32 192(%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 8), %r10
		leaq (%r10, %r8, 8), %r11

	    vpbroadcastq (%rdi), %zmm24
	    vpbroadcastq (%rdi, %r8, 8), %zmm25
	    vpbroadcastq (%r10, %r8, 8), %zmm26
	    vpbroadcastq (%r11, %r8, 8), %zmm27

        # mul
        vpmullq %zmm24, %zmm28, %zmm16
        vpmullq %zmm25, %zmm28, %zmm17
        vpmullq %zmm26, %zmm28, %zmm18
        vpmullq %zmm27, %zmm28, %zmm19

        vpmullq %zmm24, %zmm29, %zmm20
        vpmullq %zmm25, %zmm29, %zmm21
        vpmullq %zmm26, %zmm29, %zmm22
        vpmullq %zmm27, %zmm29, %zmm23

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1
		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3

		vpaddq %zmm4, %zmm20, %zmm4
		vpaddq %zmm5, %zmm21, %zmm5
		vpaddq %zmm6, %zmm22, %zmm6
		vpaddq %zmm7, %zmm23, %zmm7
        
        # mul
        vpmullq %zmm24, %zmm30, %zmm16
        vpmullq %zmm25, %zmm30, %zmm17
        vpmullq %zmm26, %zmm30, %zmm18
        vpmullq %zmm27, %zmm30, %zmm19

        vpmullq %zmm24, %zmm31, %zmm20
        vpmullq %zmm25, %zmm31, %zmm21
        vpmullq %zmm26, %zmm31, %zmm22
        vpmullq %zmm27, %zmm31, %zmm23

		# add
		vpaddq %zmm8, %zmm16, %zmm8
		vpaddq %zmm9, %zmm17, %zmm9
		vpaddq %zmm10, %zmm18, %zmm10
		vpaddq %zmm11, %zmm19, %zmm11

		vpaddq %zmm12, %zmm20, %zmm12
		vpaddq %zmm13, %zmm21, %zmm13
		vpaddq %zmm14, %zmm22, %zmm14
		vpaddq %zmm15, %zmm23, %zmm15

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_4_unr_32:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_4_unr_32

	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 8)
	vmovdqu32 %zmm2, (%r10, %r9, 8)
	vmovdqu32 %zmm3, (%r11, %r9, 8)

	vmovdqu32 %zmm4, 64(%rdx)
	vmovdqu32 %zmm5, 64(%rdx, %r9, 8)
	vmovdqu32 %zmm6, 64(%r10, %r9, 8)
	vmovdqu32 %zmm7, 64(%r11, %r9, 8)

	vmovdqu32 %zmm8, 128(%rdx)
	vmovdqu32 %zmm9, 128(%rdx, %r9, 8)
	vmovdqu32 %zmm10, 128(%r10, %r9, 8)
	vmovdqu32 %zmm11, 128(%r11, %r9, 8)

	vmovdqu32 %zmm12, 192(%rdx)
	vmovdqu32 %zmm13, 192(%rdx, %r9, 8)
	vmovdqu32 %zmm14, 192(%r10, %r9, 8)
	vmovdqu32 %zmm15, 192(%r11, %r9, 8)

	ret

# 2x8 block
mini_mmm_asm_64_64_clean_2:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 8), %zmm1

	jmp .loop_cond_clean_2
	
	.loop_body_clean_2:
		# load b vecs
		vmovdqu32 (%rsi), %zmm28
		
		# load a vecs
	    vpbroadcastq (%rdi), %zmm26
	    vpbroadcastq (%rdi, %r8, 8), %zmm27

        # mul
        vpmullq %zmm26, %zmm28, %zmm16
        vpmullq %zmm27, %zmm28, %zmm17

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_2:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_2

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 8)

	ret

# 2x16 block
mini_mmm_asm_64_64_clean_2_unr_16:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 8), %zmm1

	vmovdqu32 64(%rdx), %zmm2
	vmovdqu32 64(%rdx, %r9, 8), %zmm3

	jmp .loop_cond_clean_2_unr_16
	
	.loop_body_clean_2_unr_16:
		# load b vecs
		vmovdqu32 (%rsi), %zmm28
		vmovdqu32 64(%rsi), %zmm29
		
		# load a vecs
	    vpbroadcastq (%rdi), %zmm26
	    vpbroadcastq (%rdi, %r8, 8), %zmm27

        # mul
        vpmullq %zmm26, %zmm28, %zmm16
        vpmullq %zmm27, %zmm28, %zmm17

        vpmullq %zmm26, %zmm29, %zmm18
        vpmullq %zmm27, %zmm29, %zmm19

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1

		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_2_unr_16:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_2_unr_16

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 8)

	vmovdqu32 %zmm2, 64(%rdx)
	vmovdqu32 %zmm3, 64(%rdx, %r9, 8)

	ret

# 2x32 block
mini_mmm_asm_64_64_clean_2_unr_32:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 8), %zmm1

	vmovdqu32 64(%rdx), %zmm2
	vmovdqu32 64(%rdx, %r9, 8), %zmm3

	vmovdqu32 128(%rdx), %zmm4
	vmovdqu32 128(%rdx, %r9, 8), %zmm5

	vmovdqu32 192(%rdx), %zmm6
	vmovdqu32 192(%rdx, %r9, 8), %zmm7

	jmp .loop_cond_clean_2_unr_32
	
	.loop_body_clean_2_unr_32:
		# load b vecs
		vmovdqu32 (%rsi), %zmm28
		vmovdqu32 64(%rsi), %zmm29
		vmovdqu32 128(%rsi), %zmm30
		vmovdqu32 192(%rsi), %zmm31
		
		# load a vecs
	    vpbroadcastq (%rdi), %zmm26
	    vpbroadcastq (%rdi, %r8, 8), %zmm27

        # mul
        vpmullq %zmm26, %zmm28, %zmm16
        vpmullq %zmm27, %zmm28, %zmm17

        vpmullq %zmm26, %zmm29, %zmm18
        vpmullq %zmm27, %zmm29, %zmm19

        vpmullq %zmm26, %zmm30, %zmm20
        vpmullq %zmm27, %zmm30, %zmm21

        vpmullq %zmm26, %zmm31, %zmm22
        vpmullq %zmm27, %zmm31, %zmm23

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1

		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3

		vpaddq %zmm4, %zmm20, %zmm4
		vpaddq %zmm5, %zmm21, %zmm5

		vpaddq %zmm6, %zmm22, %zmm6
		vpaddq %zmm7, %zmm23, %zmm7

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_2_unr_32:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_2_unr_32

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 8)

	vmovdqu32 %zmm2, 64(%rdx)
	vmovdqu32 %zmm3, 64(%rdx, %r9, 8)

	vmovdqu32 %zmm4, 128(%rdx)
	vmovdqu32 %zmm5, 128(%rdx, %r9, 8)

	vmovdqu32 %zmm6, 192(%rdx)
	vmovdqu32 %zmm7, 192(%rdx, %r9, 8)

	ret


# 1x8 block
mini_mmm_asm_64_64_clean_1:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	

	jmp .loop_cond_clean_1
	
	.loop_body_clean_1:
		# load b vecs
		vmovdqu32 (%rsi), %zmm24
		
		# load a vecs
	    vpbroadcastq (%rdi), %zmm23

        # mul
        vpmullq %zmm23, %zmm24, %zmm16

		# add
		vpaddq %zmm0, %zmm16, %zmm0

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_1:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_1

	vmovdqu32 %zmm0, (%rdx)

	ret

# 1x16 block
mini_mmm_asm_64_64_clean_1_unr_16:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 64(%rdx), %zmm1

	jmp .loop_cond_clean_1_unr_16
	
	.loop_body_clean_1_unr_16:
		# load b vecs
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		
		# load a vecs
	    vpbroadcastq (%rdi), %zmm23

        # mul
        vpmullq %zmm23, %zmm24, %zmm16
        vpmullq %zmm23, %zmm25, %zmm17

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_1_unr_16:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_1_unr_16

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, 64(%rdx)

	ret

# 1x32 block
mini_mmm_asm_64_64_clean_1_unr_32:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 64(%rdx), %zmm1
	vmovdqu32 128(%rdx), %zmm2
	vmovdqu32 192(%rdx), %zmm3

	jmp .loop_cond_clean_1_unr_32
	
	.loop_body_clean_1_unr_32:
		# load b vecs
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		vmovdqu32 128(%rsi), %zmm26
		vmovdqu32 192(%rsi), %zmm27
		
		# load a vecs
	    vpbroadcastq (%rdi), %zmm23

        # mul
        vpmullq %zmm23, %zmm24, %zmm16
        vpmullq %zmm23, %zmm25, %zmm17
        vpmullq %zmm23, %zmm26, %zmm18
        vpmullq %zmm23, %zmm27, %zmm19

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1
		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_1_unr_32:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_1_unr_32

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, 64(%rdx)
	vmovdqu32 %zmm2, 128(%rdx)
	vmovdqu32 %zmm3, 192(%rdx)

	ret

# 1x64 block
mini_mmm_asm_64_64_clean_1_unr_64:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 64(%rdx), %zmm1
	vmovdqu32 128(%rdx), %zmm2
	vmovdqu32 192(%rdx), %zmm3
	vmovdqu32 256(%rdx), %zmm4
	vmovdqu32 320(%rdx), %zmm5
	vmovdqu32 384(%rdx), %zmm6
	vmovdqu32 448(%rdx), %zmm7

	jmp .loop_cond_clean_1_unr_64
	
	.loop_body_clean_1_unr_64:
		# load b vecs
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		vmovdqu32 128(%rsi), %zmm26
		vmovdqu32 192(%rsi), %zmm27
		vmovdqu32 256(%rsi), %zmm28
		vmovdqu32 320(%rsi), %zmm29
		vmovdqu32 384(%rsi), %zmm30
		vmovdqu32 448(%rsi), %zmm31
		
		# load a vecs
	    vpbroadcastq (%rdi), %zmm23

        # mul
        vpmullq %zmm23, %zmm24, %zmm16
        vpmullq %zmm23, %zmm25, %zmm17
        vpmullq %zmm23, %zmm26, %zmm18
        vpmullq %zmm23, %zmm27, %zmm19
        vpmullq %zmm23, %zmm28, %zmm20
        vpmullq %zmm23, %zmm29, %zmm21
        vpmullq %zmm23, %zmm30, %zmm22
        vpmullq %zmm23, %zmm31, %zmm23

		# add
		vpaddq %zmm0, %zmm16, %zmm0
		vpaddq %zmm1, %zmm17, %zmm1
		vpaddq %zmm2, %zmm18, %zmm2
		vpaddq %zmm3, %zmm19, %zmm3
		vpaddq %zmm4, %zmm20, %zmm4
		vpaddq %zmm5, %zmm21, %zmm5
		vpaddq %zmm6, %zmm22, %zmm6
		vpaddq %zmm7, %zmm23, %zmm7

		# adjust loop vars
		decq %rcx
		addq $8, %rdi
		leaq (%rsi, %r9, 8), %rsi
	.loop_cond_clean_1_unr_64:
		# check iterations
		cmpl $0, %ecx
    jg .loop_body_clean_1_unr_64

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, 64(%rdx)
	vmovdqu32 %zmm2, 128(%rdx)
	vmovdqu32 %zmm3, 192(%rdx)
	vmovdqu32 %zmm4, 256(%rdx)
	vmovdqu32 %zmm5, 320(%rdx)
	vmovdqu32 %zmm6, 384(%rdx)
	vmovdqu32 %zmm7, 448(%rdx)

	ret
