%include "sseutils64.nasm"
 
section .data
alpha_phi dq -57.8
alpha_psi dq -47.0
beta_phi  dq -119.0
beta_psi  dq 113.0
half      dq 0.5

 
section .bss            
 
alignb 32
e       resq        1
 
section .text          
 
extern get_block
extern free_block
 
%macro  getmem  2
    mov rdi, %1
    mov rsi, %2
    call    get_block
%endmacro
 
%macro  fremem  1
    mov rdi, %1
    call    free_block
%endmacro
 
 
global rama_energy_assembly
 
rama_energy_assembly:
    push    rbp            
    mov     rbp, rsp      
    pushaq                
 
    mov     rdi, [rbp + 16]
    mov     rsi, [rbp + 24]
   
    vmovsd  xmm0, qword [alpha_phi]
    vbroadcastsd ymm0, xmm0
    vmovsd  xmm1, qword [alpha_psi]
    vbroadcastsd ymm1, xmm1
    vmovsd  xmm2, qword [beta_phi]
    vbroadcastsd ymm2, xmm2
    vmovsd  xmm3, qword [beta_psi]
    vbroadcastsd ymm3, xmm3
 
    vxorpd  ymm4, ymm4, ymm4
 
    mov     rcx, 256      
    xor     r8, r8        
 
.loop:
 
    vmovsd  xmm5, qword [rdi + r8 * 8]
    vbroadcastsd ymm5, xmm5
    vmovsd  xmm6, qword [rsi + r8 * 8]
    vbroadcastsd ymm6, xmm6
 
    vsubpd  ymm7, ymm5, ymm0
    vmulpd  ymm7, ymm7, ymm7
 
    vsubpd  ymm8, ymm6, ymm1
    vmulpd  ymm8, ymm8, ymm8
 
    vaddpd  ymm7, ymm7, ymm8
 
    vsqrtpd ymm7, ymm7
 
    vsubpd  ymm8, ymm5, ymm2
    vmulpd  ymm8, ymm8, ymm8
 
    vsubpd  ymm9, ymm6, ymm3
    vmulpd  ymm9, ymm9, ymm9
 
    vaddpd  ymm8, ymm8, ymm9
 
    vsqrtpd ymm8, ymm8
 
    vminpd  ymm7, ymm7, ymm8
 
    vmulpd  ymm7, ymm7, [half]
    vaddpd  ymm4, ymm4, ymm7
 
    inc     r8
    loop    .loop
    vmovsd  qword [rdx], xmm4
 
    popaq                
    mov     rsp, rbp      
    pop     rbp            
    ret                    
