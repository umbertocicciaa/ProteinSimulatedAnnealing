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

section .data			; Sezione contenente dati inizializzati
    

section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	e		resd		1

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

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

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global sub

V 	equ		8
W   equ     12    
res equ     16

sub:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione
		
		MOV EAX, [EBP+V]	
        MOV EBX, [EBP+W]
		MOV ECX, [EBP+res]

		MOVAPS XMM0, [EAX]
        MOVAPS XMM1, [EBX]
        SUBPS XMM1, XMM0
        MOVAPS [ECX], XMM1
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante

global sum
V 	equ		8
W   equ     12    
res equ     16


    sum:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione
		
		
		MOV EAX, [EBP+V]	
        MOV EBX, [EBP+W]
		MOV ECX, [EBP+res]

		MOVAPS XMM0, [EAX]
        MOVAPS XMM1, [EBX]
        ADDPS XMM1, XMM0
        MOVAPS [ECX], XMM1
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante


global euclidean_dist_sse

A equ 8
B equ 12
C equ 16

msg	db	'e:',32,0


	euclidean_dist_sse:

		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione


		MOV EAX, [EBP+A]
		MOV EBX, [EBP+B]
		MOV ECX, [EBP+C]
		
		; Carica i valori di v1 (x1, y1, z1) in un registro SSE
		MOVAPS   XMM0, [EAX]           ; xmm0 = v1 (x1, y1, z1, 0)
		
		; Carica i valori di v2 (x2, y2, z2) in un altro registro SSE
		MOVAPS   XMM1, [EBX]           ; xmm1 = v2 (x2, y2, z2, 0)
		
		; Calcola le differenze (x2 - x1), (y2 - y1), (z2 - z1)
		SUBPS XMM1, XMM0          ; xmm1 = v2 - v1 
		
		; Calcola il quadrato delle differenze
		MULPS XMM1, XMM1         ; xmm1 = (x2-x1)^2, (y2-y1)^2, (z2-z1)^2, 0
		MOVSS [e], XMM1
		printss e
		
		; Somma le differenze quadrate: (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2
		HADDPS XMM1, XMM1          ; xmm1 = ((x2-x1)^2 + (y2-y1)^2), (z2-z1)^2, 0, 0
		HADDPS XMM1, XMM1          ; xmm1 = ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2), 0, 0, 0
		MOVSS [e], XMM1
		printss e

		; Calcola la radice quadrata
		SQRTPS XMM1,XMM1          ; xmm1 = sqrt(xmm1)
		MOVSS [e], XMM1
		printss e
		
		; Salva il risultato (la distanza)
		MOVSS [ECX], XMM1
		MOVSS [e], XMM1
		printss e

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------


		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante
