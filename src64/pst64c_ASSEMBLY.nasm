; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli, F. Fassetti, S. Nisticò
; 12/11/2024
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
;     nasm -f elf64 pst64.nasm
;

%include "sseutils64.nasm"


section .data
alpha_phi dq -57.8
alpha_psi dq -47.0
beta_phi  dq -119.0
beta_psi  dq 113.0
half      dq 0.5
section .bss			; Sezione contenente dati non inizializzati

alignb 32
e		resq		1

section .text			; Sezione contenente il codice macchina

;

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro


global rama_energy_assembly

rama_energy_assembly:
    ; Prologo
    push    rbp             ; salva il Base Pointer
    mov     rbp, rsp        ; il Base Pointer punta al Record di Attivazione corrente
    pushaq                  ; salva i registri generali

    ; Carica i parametri
    mov     rdi, [rbp + 16] ; phi
    mov     rsi, [rbp + 24] ; psi
   
    ; Inizializza i valori costanti
    vmovsd  xmm0, qword [alpha_phi]
    vbroadcastsd ymm0, xmm0
    vmovsd  xmm1, qword [alpha_psi]
    vbroadcastsd ymm1, xmm1
    vmovsd  xmm2, qword [beta_phi]
    vbroadcastsd ymm2, xmm2
    vmovsd  xmm3, qword [beta_psi]
    vbroadcastsd ymm3, xmm3

    ; Inizializza energia a zero
    vxorpd  ymm4, ymm4, ymm4

    ; Loop per calcolare l'energia
    mov     rcx, 256        ; n = 256
    xor     r8, r8          ; i = 0

.loop:
    ; Carica phi[i] e psi[i]
    vmovsd  xmm5, qword [rdi + r8 * 8]
    vbroadcastsd ymm5, xmm5 ;phi[i]
    vmovsd  xmm6, qword [rsi + r8 * 8]
    vbroadcastsd ymm6, xmm6 ;psi[i]

    ; Calcola (phi[i] - alpha_phi)^2
    vsubpd  ymm7, ymm5, ymm0
    vmulpd  ymm7, ymm7, ymm7

    ; Calcola (psi[i] - alpha_psi)^2
    vsubpd  ymm8, ymm6, ymm1
    vmulpd  ymm8, ymm8, ymm8

    ; Somma i quadrati
    vaddpd  ymm7, ymm7, ymm8

    ; Calcola la radice quadrata
    vsqrtpd ymm7, ymm7

    ; Calcola (phi[i] - beta_phi)^2
    vsubpd  ymm8, ymm5, ymm2
    vmulpd  ymm8, ymm8, ymm8

    ; Calcola (psi[i] - beta_psi)^2
    vsubpd  ymm9, ymm6, ymm3
    vmulpd  ymm9, ymm9, ymm9

    ; Somma i quadrati
    vaddpd  ymm8, ymm8, ymm9

    ; Calcola la radice quadrata
    vsqrtpd ymm8, ymm8

    ; Trova il minimo tra alpha_dist e beta_dist
    vminpd  ymm7, ymm7, ymm8

    ; Somma 0.5 * min
    vmulpd  ymm7, ymm7, [half]
    vaddpd  ymm4, ymm4, ymm7

    ; Incrementa l'indice
    inc     r8
    loop    .loop

    
    ; Salva il risultato in rama
    vmovsd  qword [rdx], xmm4

    ; Epilogo
    popaq                   ; ripristina i registri generali
    mov     rsp, rbp        ; ripristina lo Stack Pointer
    pop     rbp             ; ripristina il Base Pointer
    ret                     ; torna alla funzione C chiamante


