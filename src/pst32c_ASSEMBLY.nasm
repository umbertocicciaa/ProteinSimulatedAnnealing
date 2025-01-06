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
    pi dq 3.14159265358979323846
    
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

global rotation

rotation:
    ; Prologo della funzione
    push    ebp
    mov     ebp, esp
    sub     esp, 16
    push    ebx
    push    esi
    push    edi
 
    ; Carica i parametri
    mov     eax, [ebp + 8]    ; axis
    mov     ebx, [ebp + 12]   ; theta
    mov     ecx, [ebp + 16]   ; rotation_matrix
 
    ; Calcola il prodotto scalare di axis con se stesso
    movaps  xmm0, [eax]
    mulps   xmm0, xmm0
    haddps  xmm0, xmm0
    haddps  xmm0, xmm0
    movss   [esp], xmm0
 
    ; Dividi axis per il prodotto scalare
    movss   xmm1, [esp]
    movaps  xmm2, [eax]
    divps   xmm2, xmm1
    movaps  [eax], xmm2
 
    ; Moltiplica axis per approx_sin(theta / 2.0f)
    movss   xmm3, [ebx]
    mulss   xmm3, dword [pi]
    divss   xmm3, dword [2.0]
    call    approx_sin
    movaps  xmm4, [eax]
    mulps   xmm2, xmm4
    movaps  [eax], xmm2
 
    ; Moltiplica axis per -1
    movaps  xmm5, [eax]
    mulps   xmm5, dword [-1.0]
    movaps  [eax], xmm5
 
    ; Calcola approx_cos(theta / 2.0f)
    movss   xmm6, [ebx]
    mulss   xmm6, dword [pi]
    divss   xmm6, dword [2.0]
    call    approx_cos
    movss   xmm7, [eax]
 
    ; Calcola i termini della matrice di rotazione
    ; rotation_matrix[0] = a * a + b * b - c * c - d * d
    movss   xmm0, xmm7
    mulss   xmm0, xmm7
    movss   xmm1, [eax]
    mulss   xmm1, [eax]
    subss   xmm0, xmm1
    movss   [ecx], xmm0
 
    ; rotation_matrix[1] = 2 * (b * c + a * d)
    movss   xmm2, [eax + 4]
    mulss   xmm2, [eax + 8]
    addss   xmm2, xmm7
    mulss   xmm2, dword [2.0]
    movss   [ecx + 4], xmm2
 
    ; rotation_matrix[2] = 2 * (b * d - a * c)
    movss   xmm3, [eax + 4]
    mulss   xmm3, [eax + 12]
    subss   xmm3, xmm7
    mulss   xmm3, dword [2.0]
    movss   [ecx + 8], xmm3
 
    ; rotation_matrix[3] = 2 * (b * c - a * d)
    movss   xmm4, [eax + 4]
    mulss   xmm4, [eax + 8]
    subss   xmm4, xmm7
    mulss   xmm4, dword [2.0]
    movss   [ecx + 12], xmm4
 
    ; rotation_matrix[4] = a * a + c * c - b * b - d * d
    movss   xmm5, xmm7
    mulss   xmm5, xmm7
    movss   xmm6, [eax + 8]
    mulss   xmm6, [eax + 8]
    subss   xmm5, xmm6
    movss   [ecx + 16], xmm5
 
    ; rotation_matrix[5] = 2 * (c * d + a * b)
    movss   xmm7, [eax + 8]
    mulss   xmm7, [eax + 12]
    addss   xmm7, xmm7
    mulss   xmm7, dword [2.0]
    movss   [ecx + 20], xmm7
 
    ; rotation_matrix[6] = 2 * (b * d + a * c)
    movss   xmm0, [eax + 4]
    mulss   xmm0, [eax + 12]
    addss   xmm0, xmm7
    mulss   xmm0, dword [2.0]
    movss   [ecx + 24], xmm0
 
    ; rotation_matrix[7] = 2 * (c * d - a * b)
    movss   xmm1, [eax + 8]
    mulss   xmm1, [eax + 12]
    subss   xmm1, xmm7
    mulss   xmm1, dword [2.0]
    movss   [ecx + 28], xmm1
 
    ; rotation_matrix[8] = a * a + d * d - b * b - c * c
    movss   xmm2, xmm7
    mulss   xmm2, xmm7
    movss   xmm3, [eax + 12]
    mulss   xmm3, [eax + 12]
    subss   xmm2, xmm3
    movss   [ecx + 32], xmm2
 
    ; Epilogo della funzione
    pop     edi
    pop     esi
    pop     ebx
    mov     esp, ebp
    pop     ebp
    ret