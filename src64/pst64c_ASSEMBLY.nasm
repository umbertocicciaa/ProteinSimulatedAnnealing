%include "sseutils64.nasm"
 
section .data
alpha_phi dq -57.8
alpha_psi dq -47.0
beta_phi  dq -119.0
beta_psi  dq 113.0
half      dq 0.5

ten       dq 10.0
four      dq 4.0

charge dq 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.5, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, -1.0
 
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
 
 global electrostatic_energy_assembly

 electrostatic_energy_assembly:
;prologo
    push    rbp            
    mov     rbp, rsp      
    pushaq      

;elaborazione
 ; Load parameters
    mov     rdi, [rbp + 16]  ; char *s
    mov     rsi, [rbp + 24]  ; int n
    mov     rdx, [rbp + 32]  ; MATRIX coords
    mov     r8, 0            ; i = 0
    vxorpd  ymm0, ymm0, ymm0 ; energy = 0

.loop_i:
    cmp     r8, rsi
    jge     .end_loop_i

    ; Load coords_c_alpha_i
    mov     r9, r8
    imul    r9, r9, 72       ; idx_i = i * 3 * 3 * 8
    vmovupd ymm1, [rdx + r9 + 24]
    

    ; Load charge_i
    movzx   eax, byte [rdi + r8]
    sub     eax, 65
    vmovsd  xmm2, qword [charge + rax * 8]
    vbroadcastsd ymm2, xmm2

    ; Inner loop
    mov     r10, r8
    inc     r10

.loop_j:
    cmp     r10, rsi
    jge     .end_loop_j

    ; Load coords_c_alpha_j
    mov     r11, r10
    imul    r11, r11, 72      ; idx_j = j * 3 * 3 * 8
    vmovupd ymm3, [rdx + r11 + 24]

    ; Load charge_j
    movzx   eax, byte [rdi + r10]
    sub     eax, 65
    vmovsd  xmm4, qword [charge + rax * 8]
    vbroadcastsd ymm4, xmm4

    ; Calculate distance
    vsubpd  ymm5, ymm3, ymm1
    vmulpd  ymm5, ymm5, ymm5
    vhaddpd ymm5, ymm5, ymm5
    vhaddpd ymm5, ymm5, ymm5
    vsqrtpd ymm5, ymm5

    ; Calculate energy contribution
    vcmppd  ymm6, ymm5, [ten], 1 ; Compare dist < 10.0
    vmulpd  ymm7, ymm2, ymm4
    vdivpd  ymm7, ymm7, ymm5
    vdivpd  ymm7, ymm7, [four]
    vandpd  ymm7, ymm7, ymm6
    vaddpd  ymm0, ymm0, ymm7

    inc     r10
    jmp     .loop_j

.end_loop_j:
    inc     r8
    jmp     .loop_i

.end_loop_i:
    ; Store result
    vmovsd  qword [rbp + 40], xmm0
 ;epilogo
    popaq                
    mov     rsp, rbp      
    pop     rbp            
    ret                    