# file mini_mmm_asm_16_32.s
# author Jonas Kühne (jonas.kuehne@proton.me)
# Contains kernels for 16->32

# args:         rdi, rsi, rdx, rcx, r8, r9
# return:       rax
# caller saved: rax, rdi, rsi, rdx, rcx, r8, r9, r10, and r11;
# callee saved: rbx, rsp, rbp, r12, r13, r14, and r15;
.text
.globl mini_mmm_asm_16_32

.globl mini_mmm_asm_16_32_clean_8
.globl mini_mmm_asm_16_32_clean_8_unr_32

.globl mini_mmm_asm_16_32_clean_4
.globl mini_mmm_asm_16_32_clean_4_unr_32
.globl mini_mmm_asm_16_32_clean_4_unr_64

.globl mini_mmm_asm_16_32_clean_2
.globl mini_mmm_asm_16_32_clean_2_unr_32
.globl mini_mmm_asm_16_32_clean_2_unr_64
.globl mini_mmm_asm_16_32_clean_2_unr_128

.globl mini_mmm_asm_16_32_clean_1
.globl mini_mmm_asm_16_32_clean_1_unr_32
.globl mini_mmm_asm_16_32_clean_1_unr_64
.globl mini_mmm_asm_16_32_clean_1_unr_128

# %rdi: a_buf + i1*P_U
# %rsi: b_buf + (j1 << 1)
# %rdx: c_buf + i1*M_U + j1
# %rcx: k_inc
# %r8: P_U
# %r9: M_U

# scratch: %r10, %r11, %rax
# 16x16 block
mini_mmm_asm_16_32:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1
	vmovdqu32 (%rdx, %r9, 8), %zmm2
	
	vmovdqu32 (%r10, %r9, 4), %zmm3
	vmovdqu32 (%r10, %r9, 8), %zmm4
	
	leaq (%rax, %r9, 8), %r10

	vmovdqu32 (%r11, %r9, 4), %zmm5
	vmovdqu32 (%r11, %r9, 8), %zmm6

	leaq (%r10, %r9, 8), %r11

	vmovdqu32 (%rax, %r9, 4), %zmm7
	vmovdqu32 (%rax, %r9, 8), %zmm8

	leaq (%r11, %r9, 8), %rax

	vmovdqu32 (%r10, %r9, 4), %zmm9
	vmovdqu32 (%r10, %r9, 8), %zmm10
	
	leaq (%rax, %r9, 8), %r10

	vmovdqu32 (%r11, %r9, 4), %zmm11
	vmovdqu32 (%r11, %r9, 8), %zmm12

	vmovdqu32 (%rax, %r9, 4), %zmm13
	vmovdqu32 (%rax, %r9, 8), %zmm14

	vmovdqu32 (%r10, %r9, 4), %zmm15

	jmp .loop_cond
	
	.loop_body:
		# load b vec
		vmovdqu32 (%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 4), %r10
		leaq (%r10, %r8, 4), %r11
		leaq (%r11, %r8, 4), %rax

		vpbroadcastd (%rdi), %zmm16
		vpbroadcastd (%rdi, %r8, 2), %zmm17
		vpbroadcastd (%rdi, %r8, 4), %zmm18
		
		vpbroadcastd (%r10, %r8, 2), %zmm19
		vpbroadcastd (%r10, %r8, 4), %zmm20
		
		leaq (%rax, %r8, 4), %r10

		vpbroadcastd (%r11, %r8, 2), %zmm21
		vpbroadcastd (%r11, %r8, 4), %zmm22

		leaq (%r10, %r8, 4), %r11

		vpbroadcastd (%rax, %r8, 2), %zmm23
		vpbroadcastd (%rax, %r8, 4), %zmm24

		leaq (%r11, %r8, 4), %rax

		vpbroadcastd (%r10, %r8, 2), %zmm25
		vpbroadcastd (%r10, %r8, 4), %zmm26
		
		leaq (%rax, %r8, 4), %r10

		vpbroadcastd (%r11, %r8, 2), %zmm27
		vpbroadcastd (%r11, %r8, 4), %zmm28

		vpbroadcastd (%rax, %r8, 2), %zmm29
		vpbroadcastd (%rax, %r8, 4), %zmm30
		
		# fma
		vpdpwssd %zmm16, %zmm31, %zmm0

		vpbroadcastd (%r10, %r8, 2), %zmm16

		vpdpwssd %zmm17, %zmm31, %zmm1
		vpdpwssd %zmm18, %zmm31, %zmm2
		vpdpwssd %zmm19, %zmm31, %zmm3
		vpdpwssd %zmm20, %zmm31, %zmm4
		vpdpwssd %zmm21, %zmm31, %zmm5
		vpdpwssd %zmm22, %zmm31, %zmm6
		vpdpwssd %zmm23, %zmm31, %zmm7
		vpdpwssd %zmm24, %zmm31, %zmm8
		vpdpwssd %zmm25, %zmm31, %zmm9
		vpdpwssd %zmm26, %zmm31, %zmm10
		vpdpwssd %zmm27, %zmm31, %zmm11
		vpdpwssd %zmm28, %zmm31, %zmm12
		vpdpwssd %zmm29, %zmm31, %zmm13
		vpdpwssd %zmm30, %zmm31, %zmm14
		vpdpwssd %zmm16, %zmm31, %zmm15

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body
	
	# store c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)
	vmovdqu32 %zmm2, (%rdx, %r9, 8)
	
	vmovdqu32 %zmm3, (%r10, %r9, 4)
	vmovdqu32 %zmm4, (%r10, %r9, 8)
	
	leaq (%rax, %r9, 8), %r10

	vmovdqu32 %zmm5, (%r11, %r9, 4)
	vmovdqu32 %zmm6, (%r11, %r9, 8)

	leaq (%r10, %r9, 8), %r11

	vmovdqu32 %zmm7, (%rax, %r9, 4)
	vmovdqu32 %zmm8, (%rax, %r9, 8)

	leaq (%r11, %r9, 8), %rax

	vmovdqu32 %zmm9, (%r10, %r9, 4)
	vmovdqu32 %zmm10, (%r10, %r9, 8)
	
	leaq (%rax, %r9, 8), %r10

	vmovdqu32 %zmm11, (%r11, %r9, 4)
	vmovdqu32 %zmm12, (%r11, %r9, 8)

	vmovdqu32 %zmm13, (%rax, %r9, 4)
	vmovdqu32 %zmm14, (%rax, %r9, 8)

	vmovdqu32 %zmm15, (%r10, %r9, 4)

	ret

# 8x16 block
mini_mmm_asm_16_32_clean_8:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1
	vmovdqu32 (%rdx, %r9, 8), %zmm2
	vmovdqu32 (%r10, %r9, 4), %zmm3
	vmovdqu32 (%r10, %r9, 8), %zmm4
	vmovdqu32 (%r11, %r9, 4), %zmm5
	vmovdqu32 (%r11, %r9, 8), %zmm6
	vmovdqu32 (%rax, %r9, 4), %zmm7

	jmp .loop_cond_clean_8
	
	.loop_body_clean_8:
		# load b vec
		vmovdqu32 (%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 4), %r10
		leaq (%r10, %r8, 4), %r11
		leaq (%r11, %r8, 4), %rax

		vpbroadcastd (%rdi), %zmm16
		vpbroadcastd (%rdi, %r8, 2), %zmm17
		vpbroadcastd (%rdi, %r8, 4), %zmm18
		vpbroadcastd (%r10, %r8, 2), %zmm19
		vpbroadcastd (%r10, %r8, 4), %zmm20
		vpbroadcastd (%r11, %r8, 2), %zmm21
		vpbroadcastd (%r11, %r8, 4), %zmm22
		vpbroadcastd (%rax, %r8, 2), %zmm23
		
		# fma
		vpdpwssd %zmm16, %zmm31, %zmm0
		vpdpwssd %zmm17, %zmm31, %zmm1
		vpdpwssd %zmm18, %zmm31, %zmm2
		vpdpwssd %zmm19, %zmm31, %zmm3
		vpdpwssd %zmm20, %zmm31, %zmm4
		vpdpwssd %zmm21, %zmm31, %zmm5
		vpdpwssd %zmm22, %zmm31, %zmm6
		vpdpwssd %zmm23, %zmm31, %zmm7

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_8:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_8
	
	# store c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)
	vmovdqu32 %zmm2, (%rdx, %r9, 8)
	vmovdqu32 %zmm3, (%r10, %r9, 4)
	vmovdqu32 %zmm4, (%r10, %r9, 8)
	vmovdqu32 %zmm5, (%r11, %r9, 4)
	vmovdqu32 %zmm6, (%r11, %r9, 8)
	vmovdqu32 %zmm7, (%rax, %r9, 4)

	ret

# 8x32 block
mini_mmm_asm_16_32_clean_8_unr_32:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1
	vmovdqu32 (%rdx, %r9, 8), %zmm2
	vmovdqu32 (%r10, %r9, 4), %zmm3
	vmovdqu32 (%r10, %r9, 8), %zmm4
	vmovdqu32 (%r11, %r9, 4), %zmm5
	vmovdqu32 (%r11, %r9, 8), %zmm6
	vmovdqu32 (%rax, %r9, 4), %zmm7

	vmovdqu32 64(%rdx), %zmm8	
	vmovdqu32 64(%rdx, %r9, 4), %zmm9
	vmovdqu32 64(%rdx, %r9, 8), %zmm10
	vmovdqu32 64(%r10, %r9, 4), %zmm11
	vmovdqu32 64(%r10, %r9, 8), %zmm12
	vmovdqu32 64(%r11, %r9, 4), %zmm13
	vmovdqu32 64(%r11, %r9, 8), %zmm14
	vmovdqu32 64(%rax, %r9, 4), %zmm15

	jmp .loop_cond_clean_8_unr_32
	
	.loop_body_clean_8_unr_32:
		# load b vec
		vmovdqu32 (%rsi), %zmm30
		vmovdqu32 64(%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 4), %r10
		leaq (%r10, %r8, 4), %r11
		leaq (%r11, %r8, 4), %rax

		vpbroadcastd (%rdi), %zmm16
		vpbroadcastd (%rdi, %r8, 2), %zmm17
		vpbroadcastd (%rdi, %r8, 4), %zmm18
		vpbroadcastd (%r10, %r8, 2), %zmm19
		vpbroadcastd (%r10, %r8, 4), %zmm20
		vpbroadcastd (%r11, %r8, 2), %zmm21
		vpbroadcastd (%r11, %r8, 4), %zmm22
		vpbroadcastd (%rax, %r8, 2), %zmm23
		
		# fma
		vpdpwssd %zmm16, %zmm30, %zmm0
		vpdpwssd %zmm17, %zmm30, %zmm1
		vpdpwssd %zmm18, %zmm30, %zmm2
		vpdpwssd %zmm19, %zmm30, %zmm3
		vpdpwssd %zmm20, %zmm30, %zmm4
		vpdpwssd %zmm21, %zmm30, %zmm5
		vpdpwssd %zmm22, %zmm30, %zmm6
		vpdpwssd %zmm23, %zmm30, %zmm7

		vpdpwssd %zmm16, %zmm31, %zmm8
		vpdpwssd %zmm17, %zmm31, %zmm9
		vpdpwssd %zmm18, %zmm31, %zmm10
		vpdpwssd %zmm19, %zmm31, %zmm11
		vpdpwssd %zmm20, %zmm31, %zmm12
		vpdpwssd %zmm21, %zmm31, %zmm13
		vpdpwssd %zmm22, %zmm31, %zmm14
		vpdpwssd %zmm23, %zmm31, %zmm15

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_8_unr_32:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_8_unr_32
	
	# store c vecs
	leaq (%rdx, %r9, 8), %r10
	leaq (%r10, %r9, 8), %r11
	leaq (%r11, %r9, 8), %rax

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)
	vmovdqu32 %zmm2, (%rdx, %r9, 8)
	vmovdqu32 %zmm3, (%r10, %r9, 4)
	vmovdqu32 %zmm4, (%r10, %r9, 8)
	vmovdqu32 %zmm5, (%r11, %r9, 4)
	vmovdqu32 %zmm6, (%r11, %r9, 8)
	vmovdqu32 %zmm7, (%rax, %r9, 4)

	vmovdqu32 %zmm8, 64(%rdx)
	vmovdqu32 %zmm9, 64(%rdx, %r9, 4)
	vmovdqu32 %zmm10, 64(%rdx, %r9, 8)
	vmovdqu32 %zmm11, 64(%r10, %r9, 4)
	vmovdqu32 %zmm12, 64(%r10, %r9, 8)
	vmovdqu32 %zmm13, 64(%r11, %r9, 4)
	vmovdqu32 %zmm14, 64(%r11, %r9, 8)
	vmovdqu32 %zmm15, 64(%rax, %r9, 4)

	ret

# 4x16 block
mini_mmm_asm_16_32_clean_4:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1
	vmovdqu32 (%rdx, %r9, 8), %zmm2
	vmovdqu32 (%r10, %r9, 4), %zmm3

	jmp .loop_cond_clean_4
	
	.loop_body_clean_4:
		# load b vec
		vmovdqu32 (%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 4), %r10

		vpbroadcastd (%rdi), %zmm16
		vpbroadcastd (%rdi, %r8, 2), %zmm17
		vpbroadcastd (%rdi, %r8, 4), %zmm18
		vpbroadcastd (%r10, %r8, 2), %zmm19
		
		# fma
		vpdpwssd %zmm16, %zmm31, %zmm0
		vpdpwssd %zmm17, %zmm31, %zmm1
		vpdpwssd %zmm18, %zmm31, %zmm2
		vpdpwssd %zmm19, %zmm31, %zmm3

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_4:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_4
	
	# store c vecs
	leaq (%rdx, %r9, 8), %r10

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)
	vmovdqu32 %zmm2, (%rdx, %r9, 8)
	vmovdqu32 %zmm3, (%r10, %r9, 4)

	ret


# 4x32 block
mini_mmm_asm_16_32_clean_4_unr_32:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1
	vmovdqu32 (%rdx, %r9, 8), %zmm2
	vmovdqu32 (%r10, %r9, 4), %zmm3

	vmovdqu32 64(%rdx), %zmm4
	vmovdqu32 64(%rdx, %r9, 4), %zmm5
	vmovdqu32 64(%rdx, %r9, 8), %zmm6
	vmovdqu32 64(%r10, %r9, 4), %zmm7

	jmp .loop_cond_clean_4_unr_32
	
	.loop_body_clean_4_unr_32:
		# load b vec
		vmovdqu32 (%rsi), %zmm30
		vmovdqu32 64(%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 4), %r10

		vpbroadcastd (%rdi), %zmm16
		vpbroadcastd (%rdi, %r8, 2), %zmm17
		vpbroadcastd (%rdi, %r8, 4), %zmm18
		vpbroadcastd (%r10, %r8, 2), %zmm19
		
		# fma
		vpdpwssd %zmm16, %zmm30, %zmm0
		vpdpwssd %zmm17, %zmm30, %zmm1
		vpdpwssd %zmm18, %zmm30, %zmm2
		vpdpwssd %zmm19, %zmm30, %zmm3

		vpdpwssd %zmm16, %zmm31, %zmm4
		vpdpwssd %zmm17, %zmm31, %zmm5
		vpdpwssd %zmm18, %zmm31, %zmm6
		vpdpwssd %zmm19, %zmm31, %zmm7

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_4_unr_32:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_4_unr_32
	
	# store c vecs
	leaq (%rdx, %r9, 8), %r10

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)
	vmovdqu32 %zmm2, (%rdx, %r9, 8)
	vmovdqu32 %zmm3, (%r10, %r9, 4)

	vmovdqu32 %zmm4, 64(%rdx)
	vmovdqu32 %zmm5, 64(%rdx, %r9, 4)
	vmovdqu32 %zmm6, 64(%rdx, %r9, 8)
	vmovdqu32 %zmm7, 64(%r10, %r9, 4)

	ret

# 4x64 block
mini_mmm_asm_16_32_clean_4_unr_64:
	# load c vecs
	leaq (%rdx, %r9, 8), %r10

	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1
	vmovdqu32 (%rdx, %r9, 8), %zmm2
	vmovdqu32 (%r10, %r9, 4), %zmm3

	vmovdqu32 64(%rdx), %zmm4
	vmovdqu32 64(%rdx, %r9, 4), %zmm5
	vmovdqu32 64(%rdx, %r9, 8), %zmm6
	vmovdqu32 64(%r10, %r9, 4), %zmm7

	vmovdqu32 128(%rdx), %zmm8
	vmovdqu32 128(%rdx, %r9, 4), %zmm9
	vmovdqu32 128(%rdx, %r9, 8), %zmm10
	vmovdqu32 128(%r10, %r9, 4), %zmm11

	vmovdqu32 192(%rdx), %zmm12
	vmovdqu32 192(%rdx, %r9, 4), %zmm13
	vmovdqu32 192(%rdx, %r9, 8), %zmm14
	vmovdqu32 192(%r10, %r9, 4), %zmm15

	jmp .loop_cond_clean_4_unr_64
	
	.loop_body_clean_4_unr_64:
		# load b vec
		vmovdqu32 (%rsi), %zmm28
		vmovdqu32 64(%rsi), %zmm29
		vmovdqu32 128(%rsi), %zmm30
		vmovdqu32 192(%rsi), %zmm31
		
		# load a vecs
		leaq (%rdi, %r8, 4), %r10

		vpbroadcastd (%rdi), %zmm16
		vpbroadcastd (%rdi, %r8, 2), %zmm17
		vpbroadcastd (%rdi, %r8, 4), %zmm18
		vpbroadcastd (%r10, %r8, 2), %zmm19
		
		# fma
		vpdpwssd %zmm16, %zmm28, %zmm0
		vpdpwssd %zmm17, %zmm28, %zmm1
		vpdpwssd %zmm18, %zmm28, %zmm2
		vpdpwssd %zmm19, %zmm28, %zmm3

		vpdpwssd %zmm16, %zmm29, %zmm4
		vpdpwssd %zmm17, %zmm29, %zmm5
		vpdpwssd %zmm18, %zmm29, %zmm6
		vpdpwssd %zmm19, %zmm29, %zmm7

		vpdpwssd %zmm16, %zmm30, %zmm8
		vpdpwssd %zmm17, %zmm30, %zmm9
		vpdpwssd %zmm18, %zmm30, %zmm10
		vpdpwssd %zmm19, %zmm30, %zmm11

		vpdpwssd %zmm16, %zmm31, %zmm12
		vpdpwssd %zmm17, %zmm31, %zmm13
		vpdpwssd %zmm18, %zmm31, %zmm14
		vpdpwssd %zmm19, %zmm31, %zmm15

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_4_unr_64:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_4_unr_64
	
	# store c vecs
	leaq (%rdx, %r9, 8), %r10

	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)
	vmovdqu32 %zmm2, (%rdx, %r9, 8)
	vmovdqu32 %zmm3, (%r10, %r9, 4)

	vmovdqu32 %zmm4, 64(%rdx)
	vmovdqu32 %zmm5, 64(%rdx, %r9, 4)
	vmovdqu32 %zmm6, 64(%rdx, %r9, 8)
	vmovdqu32 %zmm7, 64(%r10, %r9, 4)

	vmovdqu32 %zmm8, 128(%rdx)
	vmovdqu32 %zmm9, 128(%rdx, %r9, 4)
	vmovdqu32 %zmm10, 128(%rdx, %r9, 8)
	vmovdqu32 %zmm11, 128(%r10, %r9, 4)

	vmovdqu32 %zmm12, 192(%rdx)
	vmovdqu32 %zmm13, 192(%rdx, %r9, 4)
	vmovdqu32 %zmm14, 192(%rdx, %r9, 8)
	vmovdqu32 %zmm15, 192(%r10, %r9, 4)

	ret

# 2x16 block
mini_mmm_asm_16_32_clean_2:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1

	jmp .loop_cond_clean_2
	
	.loop_body_clean_2:
		# load b vec
		vmovdqu32 (%rsi), %zmm31
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm16
		vpbroadcastd (%rdi, %r8, 2), %zmm17
		
		# fma
		vpdpwssd %zmm16, %zmm31, %zmm0
		vpdpwssd %zmm17, %zmm31, %zmm1


		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_2:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_2
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)

	ret

# 2x32 block
mini_mmm_asm_16_32_clean_2_unr_32:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1

	vmovdqu32 64(%rdx), %zmm2
	vmovdqu32 64(%rdx, %r9, 4), %zmm3

	jmp .loop_cond_clean_2_unr_32
	
	.loop_body_clean_2_unr_32:
		# load b vec
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm20
		vpbroadcastd (%rdi, %r8, 2), %zmm21
		
		# fma
		vpdpwssd %zmm20, %zmm24, %zmm0
		vpdpwssd %zmm21, %zmm24, %zmm1

		vpdpwssd %zmm20, %zmm25, %zmm2
		vpdpwssd %zmm21, %zmm25, %zmm3

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_2_unr_32:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_2_unr_32
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)

	vmovdqu32 %zmm2, 64(%rdx)
	vmovdqu32 %zmm3, 64(%rdx, %r9, 4)

	ret

# 2x64 block
mini_mmm_asm_16_32_clean_2_unr_64:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1

	vmovdqu32 64(%rdx), %zmm2
	vmovdqu32 64(%rdx, %r9, 4), %zmm3

	vmovdqu32 128(%rdx), %zmm4
	vmovdqu32 128(%rdx, %r9, 4), %zmm5

	vmovdqu32 192(%rdx), %zmm6
	vmovdqu32 192(%rdx, %r9, 4), %zmm7

	jmp .loop_cond_clean_2_unr_64
	
	.loop_body_clean_2_unr_64:
		# load b vec
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		vmovdqu32 128(%rsi), %zmm26
		vmovdqu32 192(%rsi), %zmm27
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm20
		vpbroadcastd (%rdi, %r8, 2), %zmm21
		
		# fma
		vpdpwssd %zmm20, %zmm24, %zmm0
		vpdpwssd %zmm21, %zmm24, %zmm1

		vpdpwssd %zmm20, %zmm25, %zmm2
		vpdpwssd %zmm21, %zmm25, %zmm3

		vpdpwssd %zmm20, %zmm26, %zmm4
		vpdpwssd %zmm21, %zmm26, %zmm5

		vpdpwssd %zmm20, %zmm27, %zmm6
		vpdpwssd %zmm21, %zmm27, %zmm7

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_2_unr_64:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_2_unr_64
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)

	vmovdqu32 %zmm2, 64(%rdx)
	vmovdqu32 %zmm3, 64(%rdx, %r9, 4)

	vmovdqu32 %zmm4, 128(%rdx)
	vmovdqu32 %zmm5, 128(%rdx, %r9, 4)

	vmovdqu32 %zmm6, 192(%rdx)
	vmovdqu32 %zmm7, 192(%rdx, %r9, 4)

	ret

# 2x128 block
mini_mmm_asm_16_32_clean_2_unr_128:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 (%rdx, %r9, 4), %zmm1

	vmovdqu32 64(%rdx), %zmm2
	vmovdqu32 64(%rdx, %r9, 4), %zmm3

	vmovdqu32 128(%rdx), %zmm4
	vmovdqu32 128(%rdx, %r9, 4), %zmm5

	vmovdqu32 192(%rdx), %zmm6
	vmovdqu32 192(%rdx, %r9, 4), %zmm7

	vmovdqu32 256(%rdx), %zmm8
	vmovdqu32 256(%rdx, %r9, 4), %zmm9

	vmovdqu32 320(%rdx), %zmm10
	vmovdqu32 320(%rdx, %r9, 4), %zmm11

	vmovdqu32 384(%rdx), %zmm12
	vmovdqu32 384(%rdx, %r9, 4), %zmm13

	vmovdqu32 448(%rdx), %zmm14
	vmovdqu32 448(%rdx, %r9, 4), %zmm15

	jmp .loop_cond_clean_2_unr_128
	
	.loop_body_clean_2_unr_128:
		# load b vec
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		vmovdqu32 128(%rsi), %zmm26
		vmovdqu32 192(%rsi), %zmm27
		vmovdqu32 256(%rsi), %zmm28
		vmovdqu32 320(%rsi), %zmm29
		vmovdqu32 384(%rsi), %zmm30
		vmovdqu32 448(%rsi), %zmm31
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm20
		vpbroadcastd (%rdi, %r8, 2), %zmm21
		
		# fma
		vpdpwssd %zmm20, %zmm24, %zmm0
		vpdpwssd %zmm21, %zmm24, %zmm1

		vpdpwssd %zmm20, %zmm25, %zmm2
		vpdpwssd %zmm21, %zmm25, %zmm3

		vpdpwssd %zmm20, %zmm26, %zmm4
		vpdpwssd %zmm21, %zmm26, %zmm5

		vpdpwssd %zmm20, %zmm27, %zmm6
		vpdpwssd %zmm21, %zmm27, %zmm7

		vpdpwssd %zmm20, %zmm28, %zmm8
		vpdpwssd %zmm21, %zmm28, %zmm9

		vpdpwssd %zmm20, %zmm29, %zmm10
		vpdpwssd %zmm21, %zmm29, %zmm11

		vpdpwssd %zmm20, %zmm30, %zmm12
		vpdpwssd %zmm21, %zmm30, %zmm13

		vpdpwssd %zmm20, %zmm31, %zmm14
		vpdpwssd %zmm21, %zmm31, %zmm15

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_2_unr_128:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_2_unr_128
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm1, (%rdx, %r9, 4)

	vmovdqu32 %zmm2, 64(%rdx)
	vmovdqu32 %zmm3, 64(%rdx, %r9, 4)

	vmovdqu32 %zmm4, 128(%rdx)
	vmovdqu32 %zmm5, 128(%rdx, %r9, 4)

	vmovdqu32 %zmm6, 192(%rdx)
	vmovdqu32 %zmm7, 192(%rdx, %r9, 4)

	vmovdqu32 %zmm8, 256(%rdx)
	vmovdqu32 %zmm9, 256(%rdx, %r9, 4)

	vmovdqu32 %zmm10, 320(%rdx)
	vmovdqu32 %zmm11, 320(%rdx, %r9, 4)

	vmovdqu32 %zmm12, 384(%rdx)
	vmovdqu32 %zmm13, 384(%rdx, %r9, 4)

	vmovdqu32 %zmm14, 448(%rdx)
	vmovdqu32 %zmm15, 448(%rdx, %r9, 4)

	ret

# 1x16 block
mini_mmm_asm_16_32_clean_1:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	

	jmp .loop_cond_clean_1
	
	.loop_body_clean_1:
		# load b vec
		vmovdqu32 (%rsi), %zmm31
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm16
		
		# fma
		vpdpwssd %zmm16, %zmm31, %zmm0

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_1:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_1
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)

	ret


# 1x32 block
mini_mmm_asm_16_32_clean_1_unr_32:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 64(%rdx), %zmm2

	jmp .loop_cond_clean_1_unr_32
	
	.loop_body_clean_1_unr_32:
		# load b vec
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm20
		
		# fma
		vpdpwssd %zmm20, %zmm24, %zmm0
		vpdpwssd %zmm20, %zmm25, %zmm2

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_1_unr_32:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_1_unr_32
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm2, 64(%rdx)

	ret

# 1x64 block
mini_mmm_asm_16_32_clean_1_unr_64:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 64(%rdx), %zmm2
	vmovdqu32 128(%rdx), %zmm4
	vmovdqu32 192(%rdx), %zmm6

	jmp .loop_cond_clean_1_unr_64
	
	.loop_body_clean_1_unr_64:
		# load b vec
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		vmovdqu32 128(%rsi), %zmm26
		vmovdqu32 192(%rsi), %zmm27
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm20
		
		# fma
		vpdpwssd %zmm20, %zmm24, %zmm0
		vpdpwssd %zmm20, %zmm25, %zmm2
		vpdpwssd %zmm20, %zmm26, %zmm4
		vpdpwssd %zmm20, %zmm27, %zmm6

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_1_unr_64:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_1_unr_64
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm2, 64(%rdx)
	vmovdqu32 %zmm4, 128(%rdx)
	vmovdqu32 %zmm6, 192(%rdx)

	ret

# 1x128 block
mini_mmm_asm_16_32_clean_1_unr_128:
	# load c vecs
	vmovdqu32 (%rdx), %zmm0	
	vmovdqu32 64(%rdx), %zmm2
	vmovdqu32 128(%rdx), %zmm4
	vmovdqu32 192(%rdx), %zmm6
	vmovdqu32 256(%rdx), %zmm8
	vmovdqu32 320(%rdx), %zmm10
	vmovdqu32 384(%rdx), %zmm12
	vmovdqu32 448(%rdx), %zmm14

	jmp .loop_cond_clean_1_unr_128
	
	.loop_body_clean_1_unr_128:
		# load b vec
		vmovdqu32 (%rsi), %zmm24
		vmovdqu32 64(%rsi), %zmm25
		vmovdqu32 128(%rsi), %zmm26
		vmovdqu32 192(%rsi), %zmm27
		vmovdqu32 256(%rsi), %zmm28
		vmovdqu32 320(%rsi), %zmm29
		vmovdqu32 384(%rsi), %zmm30
		vmovdqu32 448(%rsi), %zmm31
		
		# load a vecs
		vpbroadcastd (%rdi), %zmm20
		
		# fma
		vpdpwssd %zmm20, %zmm24, %zmm0
		vpdpwssd %zmm20, %zmm25, %zmm2
		vpdpwssd %zmm20, %zmm26, %zmm4
		vpdpwssd %zmm20, %zmm27, %zmm6
		vpdpwssd %zmm20, %zmm28, %zmm8
		vpdpwssd %zmm20, %zmm29, %zmm10
		vpdpwssd %zmm20, %zmm30, %zmm12
		vpdpwssd %zmm20, %zmm31, %zmm14

		# adjust loop vars
		subq $2, %rcx
		addq $4, %rdi
		leaq (%rsi, %r9, 4), %rsi
	.loop_cond_clean_1_unr_128:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_clean_1_unr_128
	
	# store c vecs
	vmovdqu32 %zmm0, (%rdx)
	vmovdqu32 %zmm2, 64(%rdx)
	vmovdqu32 %zmm4, 128(%rdx)
	vmovdqu32 %zmm6, 192(%rdx)
	vmovdqu32 %zmm8, 256(%rdx)
	vmovdqu32 %zmm10, 320(%rdx)
	vmovdqu32 %zmm12, 384(%rdx)
	vmovdqu32 %zmm14, 448(%rdx)

	ret
