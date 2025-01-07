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
    
    hydrophobicity dd 1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1


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

global hydrophobic_energy_assembly

hydrophobic_energy_assembly:
    push ebp
    mov ebp, esp
    push ebx
    push esi
    push edi

    mov eax, [ebp + 8]  ; s
    mov ebx, [ebp + 12] ; n
    mov ecx, [ebp + 16] ; coords
    mov edx, [ebp + 20] ; result

    xorps xmm0, xmm0    ; energy = 0

    mov esi, 0          ; i = 0
    .outer_loop:
        cmp esi, ebx
        jge .end_outer_loop

        mov edi, esi
        add edi, 1      ; j = i + 1

        mov eax, ecx
        add eax, esi
        shl eax, 2
        add eax, ecx
        movaps xmm1, [eax + 12] ; coords_c_alpha_i

        .inner_loop:
            cmp edi, ebx
            jge .end_inner_loop

            mov eax, ecx
            add eax, edi
            shl eax, 2
            add eax, ecx
            movaps xmm2, [eax + 12] ; coords_c_alpha_j

            movaps xmm3, xmm2
            subps xmm3, xmm1
            mulps xmm3, xmm3
            movhlps xmm4, xmm3
            addps xmm3, xmm4
            movaps xmm4, xmm3
            shufps xmm4, xmm4, 1
            addss xmm3, xmm4
            sqrtss xmm3, xmm3                ; dist = sqrtf(...)

            movzx edx, byte [eax + esi]
            sub edx, 65
            movss xmm4, [hydrophobicity + edx * 4]

            movzx edx, byte [eax + edi]
            sub edx, 65
            movss xmm5, [hydrophobicity + edx * 4]

            mulss xmm4, xmm5
            divss xmm4, xmm3

            addss xmm0, xmm4

            add edi, 1
            jmp .inner_loop

        .end_inner_loop:
        add esi, 1
        jmp .outer_loop

    .end_outer_loop:
    movss [edx], xmm0

    pop edi
    pop esi
    pop ebx
    mov esp, ebp
    pop ebp
    ret