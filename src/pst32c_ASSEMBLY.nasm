; ---------------------------------------------------------
; Predizione struttura terziaria con istruzioni SSE a 32 bit
; ---------------------------------------------------------

; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data	
	;rana energy
    alpha_phi dd -57.8
    alpha_psi dd -47.0
    beta_phi dd -119.0
    beta_psi dd 113.0
    half dd 0.5
    
section .bss			
	alignb 16
	e		resd		1

section .text			


extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro


global rama_energy_assembly

rama_energy_assembly:
    push ebp
    mov ebp, esp
    push ebx
    push esi
    push edi

    mov eax, [ebp + 8]  
    mov ebx, [ebp + 12] 
    mov ecx, [ebp + 16] 

    xorps xmm0, xmm0    

    movss xmm1, [alpha_phi]
    movss xmm2, [alpha_psi]
    movss xmm3, [beta_phi]
    movss xmm4, [beta_psi]
    movss xmm5, [half]

    mov esi, 0
    .loop:
        movss xmm6, [eax + esi * 4]
        movss xmm7, [ebx + esi * 4]

        movaps xmm0, xmm6
        subss xmm0, xmm1
        mulss xmm0, xmm0

        movaps xmm2, xmm7
        subss xmm2, xmm2
        mulss xmm2, xmm2

        addss xmm0, xmm2
        sqrtss xmm0, xmm0

        movaps xmm2, xmm6
        subss xmm2, xmm3
        mulss xmm2, xmm2

        movaps xmm3, xmm7
        subss xmm3, xmm4
        mulss xmm3, xmm3

        addss xmm2, xmm3
        sqrtss xmm2, xmm2

        minss xmm0, xmm2
        mulss xmm0, xmm5
        addss xmm0, xmm0

        add esi, 1
        cmp esi, 256
        jl .loop

    movss [ecx], xmm0

    pop edi
    pop esi
    pop ebx
    mov esp, ebp
    pop ebp
    ret