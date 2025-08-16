# file func_amx.s
# author Jonas Kühne (jonas.kuehne@proton.me)
# Contains asm functions to benchmark the amx mmm instruction

# args:         rdi, rsi, rdx, rcx, r8, r9
# return:       rax
# caller saved: rax, rdi, rsi, rdx, rcx, r8, r9, r10, and r11;
# callee saved: rbx, rsp, rbp, r12, r13, r14, and r15;
.text
.globl peak_func
.globl latency_func


# %rdi: how many runs
# returns how many instructions per run
peak_func:
    # 4*6*2*16*16*64
    movq $786432, %rax
	jmp .loop_cond_peak
	.loop_body_peak:
	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm1, %tmm3
	tdpbuud %tmm0, %tmm1, %tmm4
	tdpbuud %tmm0, %tmm1, %tmm5
	tdpbuud %tmm0, %tmm1, %tmm6
	tdpbuud %tmm0, %tmm1, %tmm7

	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm1, %tmm3
	tdpbuud %tmm0, %tmm1, %tmm4
	tdpbuud %tmm0, %tmm1, %tmm5
	tdpbuud %tmm0, %tmm1, %tmm6
	tdpbuud %tmm0, %tmm1, %tmm7

	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm1, %tmm3
	tdpbuud %tmm0, %tmm1, %tmm4
	tdpbuud %tmm0, %tmm1, %tmm5
	tdpbuud %tmm0, %tmm1, %tmm6
	tdpbuud %tmm0, %tmm1, %tmm7

	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm1, %tmm3
	tdpbuud %tmm0, %tmm1, %tmm4
	tdpbuud %tmm0, %tmm1, %tmm5
	tdpbuud %tmm0, %tmm1, %tmm6
	tdpbuud %tmm0, %tmm1, %tmm7

		decq %rdi
	.loop_cond_peak:
		cmpq $0, %rdi
		jg .loop_body_peak
		
	ret

latency_func:
    # 4*7
    movq $28, %rax
	jmp .loop_cond_latency
	.loop_body_latency:
	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm2, %tmm3
	tdpbuud %tmm0, %tmm3, %tmm4
	tdpbuud %tmm0, %tmm4, %tmm5
	tdpbuud %tmm0, %tmm5, %tmm6
	tdpbuud %tmm0, %tmm6, %tmm7
	tdpbuud %tmm0, %tmm7, %tmm1

	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm2, %tmm3
	tdpbuud %tmm0, %tmm3, %tmm4
	tdpbuud %tmm0, %tmm4, %tmm5
	tdpbuud %tmm0, %tmm5, %tmm6
	tdpbuud %tmm0, %tmm6, %tmm7
	tdpbuud %tmm0, %tmm7, %tmm1

	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm2, %tmm3
	tdpbuud %tmm0, %tmm3, %tmm4
	tdpbuud %tmm0, %tmm4, %tmm5
	tdpbuud %tmm0, %tmm5, %tmm6
	tdpbuud %tmm0, %tmm6, %tmm7
	tdpbuud %tmm0, %tmm7, %tmm1

	tdpbuud %tmm0, %tmm1, %tmm2
	tdpbuud %tmm0, %tmm2, %tmm3
	tdpbuud %tmm0, %tmm3, %tmm4
	tdpbuud %tmm0, %tmm4, %tmm5
	tdpbuud %tmm0, %tmm5, %tmm6
	tdpbuud %tmm0, %tmm6, %tmm7
	tdpbuud %tmm0, %tmm7, %tmm1

	decq %rdi
	.loop_cond_latency:
		cmpq $0, %rdi
		jg .loop_body_latency
		
	ret
