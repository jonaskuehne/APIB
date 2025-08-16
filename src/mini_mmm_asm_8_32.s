# file mini_mmm_asm_8_32.s
# author Jonas Kühne (jonas.kuehne@proton.me)
# Contains kernels for 8->32

# args:         rdi, rsi, rdx, rcx, r8, r9
# return:       rax
# caller saved: rax, rdi, rsi, rdx, rcx, r8, r9, r10, and r11;
# callee saved: rbx, rsp, rbp, r12, r13, r14, and r15;
.text
.globl mini_mmm_asm_8_32
.globl mini_mmm_asm_8_32_unr

# %rdi: a_buf
# %rsi: b_buf
# %rdx: c_buf
# %rcx: k_inc
# %r8: P_U
# %r9: M_U << 2

# %r10: (M_U << 2) << 4
# scratch: %r10, %r11, %rax
mini_mmm_asm_8_32:
	# load c tile
	tileloadd (%rdx, %r9, 1), %tmm0
	
	movq %r9, %r10
	shlq $4, %r10

	jmp .loop_cond_unr
	
	.loop_body_unr:
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3		
		# load b tile
		tileloadd (%rsi, %r9, 1), %tmm4
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
        ################################
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3		
		# load b tile
		tileloadd (%rsi, %r9, 1), %tmm4
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
        ################################
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3		
		# load b tile
		tileloadd (%rsi, %r9, 1), %tmm4
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
        ################################
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3		
		# load b tile
		tileloadd (%rsi, %r9, 1), %tmm4
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
        ################################
	.loop_cond_unr:
		# check iterations
		cmpl $255, %ecx
        jg .loop_body_unr
	
	jmp .loop_cond_clean
	
	.loop_body_clean:
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3		
		# load b tile
		tileloadd (%rsi, %r9, 1), %tmm4
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
	.loop_cond_clean:
		# check iterations
		cmpl $0, %ecx
        jg .loop_body_clean
	# store c tile
	tilestored %tmm0, (%rdx, %r9, 1)

	ret


mini_mmm_asm_8_32_unr:
	# load c tiles
	tileloadd (%rdx, %r9, 1), %tmm0
	tileloadd 64(%rdx, %r9, 1), %tmm1
	
	movq %r9, %r10
	shlq $4, %r10

	jmp .loop_cond_unr_unr
	
	.loop_body_unr_unr:
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3
		# load b tiles
		tileloadd (%rsi, %r9, 1), %tmm4
		tileloadd 64(%rsi, %r9, 1), %tmm5
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		tdpbuud %tmm5, %tmm3, %tmm1
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
		addq %r10, %rax
        #################################
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3
		# load b tiles
		tileloadd (%rsi, %r9, 1), %tmm4
		tileloadd 64(%rsi, %r9, 1), %tmm5
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		tdpbuud %tmm5, %tmm3, %tmm1
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
		addq %r10, %rax
        #################################
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3
		# load b tiles
		tileloadd (%rsi, %r9, 1), %tmm4
		tileloadd 64(%rsi, %r9, 1), %tmm5
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		tdpbuud %tmm5, %tmm3, %tmm1
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
		addq %r10, %rax
        #################################
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3
		# load b tiles
		tileloadd (%rsi, %r9, 1), %tmm4
		tileloadd 64(%rsi, %r9, 1), %tmm5
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		tdpbuud %tmm5, %tmm3, %tmm1

		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
		addq %r10, %rax

	.loop_cond_unr_unr:
		# check iterations
		cmpl $255, %ecx
		jg .loop_body_unr_unr

	jmp .loop_cond_unr_clean
	
	.loop_body_unr_clean:
		# load a tile
		tileloadd (%rdi, %r8, 1), %tmm3
		# load b tiles
		tileloadd (%rsi, %r9, 1), %tmm4
		tileloadd 64(%rsi, %r9, 1), %tmm5
		# multiply
		tdpbuud %tmm4, %tmm3, %tmm0
		tdpbuud %tmm5, %tmm3, %tmm1
		
		# adjust indices
		subl $64, %ecx
		addq $64, %rdi
		addq %r10, %rsi
		addq %r10, %rax

	.loop_cond_unr_clean:
		# check iterations
		cmpl $0, %ecx
		jg .loop_body_unr_clean
	
	# store c tile
	tilestored %tmm0, (%rdx, %r9, 1)
	tilestored %tmm1, 64(%rdx, %r9, 1)

	ret	
	
	


