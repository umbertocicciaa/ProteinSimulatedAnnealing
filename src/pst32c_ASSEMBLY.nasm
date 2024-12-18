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

global sum
res 	equ	16
W       equ     12    
V       equ     8

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
        MOVAPS XMM0, [EAX]
        MOV EBX, [EBP+W]
        MOVAPS XMM1, [EBX]
        MOV ECX, [EBP+res]
        ADDPS XMM1, XMM0
        MOVUPS [ECX], XMM1
		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante
                                                                                      
global sub
res     equ   16 
W       equ   12
V       equ   8


    sub:
                ; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione
                ; ------------------------------------------------------------
               push            ebp             ; salva il Base Pointer
               mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
               push            ebx             ; salva i registri da preservare
               push            esi
               push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione


	MOV EAX, [EBP+V]
	MOVAPS XMM0, [EAX]
  	MOV EBX, [EBP+W]
  	MOVAPS XMM1, [EBX]
	MOV ECX, [EBP+res]
	SUBPS XMM0, XMM1
	MOVUPS [ECX], XMM0

                ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------

              pop     edi             ; ripristina i registri da preservare
              pop     esi
              pop     ebx
              mov     esp, ebp        ; ripristina lo Stack Pointer
              pop     ebp             ; ripristina il Base Pointer
              ret                     ; torna alla funzione C chiamante



global normalize 
V equ 8
S equ 12

   normalize:
                ; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione

                ; ------------------------------------------------------------
               push            ebp             ; salva il Base Pointer
               mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
               push            ebx             ; salva i registri da preservare
               push            esi
               push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione


        MOV EAX, [EBP+V]
        MOVAPS XMM0, [EAX]
        MOV EBX, [EBP+S]
        MULPS XMM0, XMM0
	HADDPS XMM0, XMM0
	HADDPS XMM0, XMM0
        MOVSS [EBX], XMM0

                ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------

              pop     edi             ; ripristina i registri da preservare
              pop     esi
              pop     ebx
              mov     esp, ebp        ; ripristina lo Stack Pointer
              pop     ebp             ; ripristina il Base Pointer
              ret                     ; torna alla funzione C chiamante


global mul 
V equ 8
W equ 12

    mul:
                ; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione
                ; ------------------------------------------------------------
               push            ebp             ; salva il Base Pointer
               mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
               push            ebx             ; salva i registri da preservare
               push            esi
               push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione


        MOV EAX, [EBP+V]
        MOVAPS XMM0, [EAX]
        MOV EBX, [EBP+W]
        MOVAPS XMM1, [EBX]
        MULPS XMM0, XMM1
        MOVAPS [EAX], XMM0

                ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------

              pop     edi             ; ripristina i registri da preservare
              pop     esi
              pop     ebx
              mov     esp, ebp        ; ripristina lo Stack Pointer
              pop     ebp             ; ripristina il Base Pointer
              ret                     ; torna alla funzione C chiamante

global euclidean_dist_sse

A equ 8
B equ 12
C equ 16

msg2	db	'e:',32,0


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
		
		; Somma le differenze quadrate: (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2
		HADDPS XMM1, XMM1          ; xmm1 = ((x2-x1)^2 + (y2-y1)^2), (z2-z1)^2, 0, 0
		HADDPS XMM1, XMM1          ; xmm1 = ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2), 0, 0, 0
		MOVSS [e], XMM1

		; Calcola la radice quadrata
		SQRTPS XMM1,XMM1          ; xmm1 = sqrt(xmm1)
		MOVSS [e], XMM1
		
		; Salva il risultato (la distanza)
		MOVSS [ECX], XMM1
		MOVSS [e], XMM1

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------


		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante
 
global div_scalare
V equ 8
W equ 12

    div_scalare:
                ; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione
                ; ------------------------------------------------------------
               push            ebp             ; salva il Base Pointer
               mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
               push            ebx             ; salva i registri da preservare
               push            esi
               push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione


        MOV EAX, [EBP+V]
        MOVAPS XMM0, [EAX]
        MOV EBX, [EBP+W]
        MOVAPS XMM1, [EBX]
        DIVPS XMM0, XMM1
        MOVAPS [EAX], XMM0

                ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------

              pop     edi             ; ripristina i registri da preservare
              pop     esi
              pop     ebx
              mov     esp, ebp        ; ripristina lo Stack Pointer
              pop     ebp             ; ripristina il Base Pointer
              ret                     ; torna alla funzione C chiamante

global alpha_beta_dist

V equ 8
W equ 12
res1 equ 16
res2 equ 20

	alpha_beta_dist:
	; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione
                ; ------------------------------------------------------------
                push            ebp             ; salva il Base Pointer
                mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
                push            ebx             ; salva i registri da preservare
                push            esi
                push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione
                MOV EAX, [EBP+V]
                MOV EBX, [EBP+W]
                MOVAPS XMM0, [EAX]
		MOVAPS XMM1, [EBX]
		MOV ECX, [EBP+res1]
		MOV EDX, [EBP+res2]
		SUBPS XMM0, XMM1
		MULPS XMM0, XMM0
                HADDPS XMM0, XMM0
                SQRTPS XMM0, XMM0
                MOVSS [ECX], XMM0
		MOVSHDUP XMM0, XMM0 	
		MOVSS [EDX], XMM0
		
  ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------


                pop     edi             ; ripristina i registri da preservare
                pop     esi
                pop     ebx
                mov     esp, ebp        ; ripristina lo Stack Pointer
                pop     ebp             ; ripristina il Base Pointer
                ret                     ; torna alla funzione C chiamante

global prodotto_scalare
V equ 8
W equ 12
RES equ 16
   
    prodotto_scalare:
                ; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione
                ; ------------------------------------------------------------
               push            ebp             ; salva il Base Pointer
               mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
               push            ebx             ; salva i registri da preservare
               push            esi
               push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione


        MOV EAX, [EBP+V]
        MOVAPS XMM0, [EAX]
        MOV EBX, [EBP+W]
        MOVAPS XMM1, [EBX]
	MOV ECX, [EBP+RES]
        MULPS XMM0, XMM1
	HADDPS XMM0, XMM0
	HADDPS XMM0, XMM0
        MOVSS [ECX], XMM0

                ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------

              pop     edi             ; ripristina i registri da preservare
              pop     esi
              pop     ebx
              mov     esp, ebp        ; ripristina lo Stack Pointer
              pop     ebp             ; ripristina il Base Pointer
              ret                     ; torna alla funzione C chiamante


global prodotto_vettore_matrice_CICLO:
    A equ 8    ;v
    B equ 12   ;m
    C equ 16   ;res

    prodotto_vettore_matrice_CICLO: 
         ; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione
                ; ------------------------------------------------------------
               push            ebp             ; salva il Base Pointer
               mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
               push            ebx             ; salva i registri da preservare
               push            esi
               push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione
        MOV eax, [EBP+A] ;vettore v
        MOV ebx, [EBP+B] ;matrice m in column-major order
        MOV edx, [EBP+C] ;vettore res

    ; Inizializza res[0..3] = 0
    xorps xmm0, xmm0
    movaps [edx], xmm0          ; res[0..3] = 0

	; Carica il vettore v in xmm2 (4 float)
    movaps xmm2, [eax]          ; xmm2 = v[0..3]

    ; Ciclo per 4 iterazioni (per ogni colonna della matrice 4x4)
    mov ecx, 0                  ; ecx = indice della colonna incrementato di 4float per volta
    mov edi, 0                  ; edi = indice di posizionamento di res[i]


.loop_start:
    ; Carica la colonna della matrice nella XMM1 (usando column-major order)
    movaps xmm1, [ebx + ecx*4]  ; Carica la colonna `ecx` della matrice nella xmm1 (4 valori consecutivi)

    ; Moltiplicazione elemento per elemento tra xmm1 (colonna della matrice) e xmm2 (vettore v)
    mulps xmm1, xmm2            ; xmm1 = xmm1 * xmm2 (moltiplica ogni elemento di xmm1 per xmm2)

    ; Somma orizzontale: somma i valori in xmm1
    haddps xmm1, xmm1           ; xmm1 = [x1+x2, x3+x4, x1+x2, x3+x4]
    haddps xmm1, xmm1           ; xmm1 = [x1+x2+x3+x4, x1+x2+x3+x4, x1+x2+x3+x4, x1+x2+x3+x4]

    ; Salva il risultato in res[ecx]
    movss [edx + edi*4], xmm1   ; Copia il risultato in res[ecx] (solo il primo valore di xmm1)

    ; Incrementa l'indice della colonna
    add ecx, 4
    inc edi
    cmp ecx, 16
    jl .loop_start
  ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------

              pop     edi             ; ripristina i registri da preservare
              pop     esi
              pop     ebx
              mov     esp, ebp        ; ripristina lo Stack Pointer
              pop     ebp             ; ripristina il Base Pointer
              ret                     ; torna alla funzione C chiamante

global prodotto_vettore_matrice:
    A equ 8  ;v
    B equ 12   ;m
    C equ 16   ;res

    prodotto_vettore_matrice: 
         ; ------------------------------------------------------------
                ; Sequenza di ingresso nella funzione
                ; ------------------------------------------------------------
               push            ebp             ; salva il Base Pointer
               mov             ebp, esp        ; il Base Pointer punta al Record di Attivazione corrente
               push            ebx             ; salva i registri da preservare
               push            esi
               push            edi
                ; ------------------------------------------------------------
                ; legge i parametri dal Record di Attivazione corrente
                ; ------------------------------------------------------------

                ; elaborazione
        MOV eax, [EBP+A] ;vettore v
        MOV ebx, [EBP+B] ;matrice m in column-major order
        MOV edx, [EBP+C] ;vettore res

    ; Inizializza res[0..3] = 0
    xorps xmm0, xmm0
    movaps [edx], xmm0          ; res[0..3] = 0

	; Carica il vettore v in xmm2 (4 float)
    movaps xmm2, [eax]          ; xmm2 = v[0..3]


    ; Carica la colonna della matrice nella XMM1 (usando column-major order)
    movaps xmm1, [ebx]  ; Carica la colonna `ecx` della matrice nella xmm1 (4 valori consecutivi)
    movaps xmm3, [ebx + 16]
    movaps xmm4, [ebx + 32]
    movaps xmm5, [ebx + 48]

    ; Moltiplicazione elemento per elemento tra xmm1 (colonna della matrice) e xmm2 (vettore v)
    mulps xmm1, xmm2            ; xmm1 = xmm1 * xmm2 (moltiplica ogni elemento di xmm1 per xmm2)
    mulps xmm3, xmm2
    mulps xmm4, xmm2
    mulps xmm5, xmm2

    ; Somma orizzontale: somma i valori in xmm1
    haddps xmm1, xmm1           ; xmm1 = [x1+x2, x3+x4, x1+x2, x3+x4]
    haddps xmm1, xmm1           ; xmm1 = [x1+x2+x3+x4, x1+x2+x3+x4, x1+x2+x3+x4, x1+x2+x3+x4]

    haddps xmm3, xmm3
    haddps xmm3, xmm3

    haddps xmm4, xmm4
    haddps xmm4, xmm4

    haddps xmm5, xmm5
    haddps xmm5, xmm5

    ; Salva il risultato in res[ecx]
    movss [edx], xmm1   ; Copia il risultato in res[ecx] (solo il primo valore di xmm1)
    movss [edx + 4], xmm3
    movss [edx + 8], xmm4
    movss [edx + 12], xmm5


  ; ------------------------------------------------------------
                ; Sequenza di uscita dalla funzione
                ; ------------------------------------------------------------

              pop     edi             ; ripristina i registri da preservare
              pop     esi
              pop     ebx
              mov     esp, ebp        ; ripristina lo Stack Pointer
              pop     ebp             ; ripristina il Base Pointer
              ret                     ; torna alla funzione C chiamante
