%include "sseutils64.nasm"
 
section .data
alpha_phi dq -57.8
alpha_psi dq -47.0
beta_phi  dq -119.0
beta_psi  dq 113.0
half      dq 0.5

align 8
ten:    dq 10.0
four:   dq 4.0
charge dq 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.5, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, -1.0

hydrophobicity dq 1.8, -1.0, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1.0, -3.9, 3.8, 1.9, -3.5, -1.0, -1.6, -3.5, -4.5, -0.8, -0.7, -1.0, 4.2, -0.9, -1.0, -1.3, -1.0
 
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


    vxorpd xmm0, xmm0, xmm0    
    
    xor r8, r8                 
outer_loop:
    cmp r8, rsi                 
    jge done
    
   
    mov r9, r8
    imul r9, 9
    add r9, 3
    
    
    vmovupd ymm1, [rdx + r9*8] 
    
  
    mov r10, r8
    inc r10
    
inner_loop:
    cmp r10, rsi               
    jge outer_loop_end
    
    
    mov r11, r10
    imul r11, 9
    add r11, 3
    
    
    vmovupd ymm2, [rdx + r11*8] 
    
    
    vsubpd ymm3, ymm2, ymm1     
    vmulpd ymm3, ymm3, ymm3     
    vhaddpd ymm3, ymm3, ymm3    
    vhaddpd ymm3, ymm3, ymm3    
    
    vsqrtpd ymm3, ymm3         
    
    
    vmovsd xmm4, [rel ten]
    vucomisd xmm3, xmm4
    jae skip_update
    

    movzx r12, byte [rdi + r8]  
    sub r12, 65                 
    movzx r13, byte [rdi + r10] 
    sub r13, 65                 
    
    vmovsd xmm4, [rel charge + r12*8]
    vmovsd xmm5, [rel charge + r13*8]
    
    
    vxorpd xmm6, xmm6, xmm6
    vucomisd xmm4, xmm6
    je skip_update
    vucomisd xmm5, xmm6
    je skip_update
    
    
    vmulsd xmm4, xmm4, xmm5    
    vmulsd xmm5, xmm3, [rel four] 
    vdivsd xmm4, xmm4, xmm5     
    vaddsd xmm0, xmm0, xmm4     
    
skip_update:
    inc r10
    jmp inner_loop
    
outer_loop_end:
    inc r8
    jmp outer_loop
    
done:
    vmovsd [rcx], xmm0          
    
    popaq                       
    mov rsp, rbp
    pop rbp
    ret


global hydrophobic_energy_assembly
hydrophobic_energy_assembly:
    push    rbp
    mov     rbp, rsp
    pushaq

    vxorpd xmm0, xmm0, xmm0   

    xor r8, r8                
outer_loop_hydro:
    cmp r8, rsi               
    jge done_hydro            

    mov r9, r8
    imul r9, 9
    add r9, 3

    vmovupd ymm1, [rdx + r9*8]  

    mov r10, r8
    inc r10                  

inner_loop_hydro:
    cmp r10, rsi              
    jge outer_loop_end_hydro  

    mov r11, r10
    imul r11, 9
    add r11, 3

    vmovupd ymm2, [rdx + r11*8] 

    vsubpd ymm3, ymm2, ymm1   
    vmulpd ymm3, ymm3, ymm3   
    vhaddpd ymm3, ymm3, ymm3  
    vhaddpd ymm3, ymm3, ymm3  

    vsqrtpd ymm3, ymm3        

    vmovsd xmm4, [rel ten]
    vucomisd xmm3, xmm4
    jae skip_update_hydro     

    movzx r12, byte [rdi + r8]  
    sub r12, 65                 
    movzx r13, byte [rdi + r10] 
    sub r13, 65                 

    vmovsd xmm4, [rel hydrophobicity + r12*8]
    vmovsd xmm5, [rel hydrophobicity + r13*8]

    vmulsd xmm4, xmm4, xmm5    
    vdivsd xmm4, xmm4, xmm3    
    vaddsd xmm0, xmm0, xmm4    

skip_update_hydro:
    inc r10
    jmp inner_loop_hydro

outer_loop_end_hydro:
    inc r8
    jmp outer_loop_hydro

done_hydro:
    vmovsd [rcx], xmm0         

    popaq
    mov rsp, rbp
    pop rbp
    ret