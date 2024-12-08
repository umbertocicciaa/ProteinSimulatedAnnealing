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

global prova

input		equ		8

msg	db	'e:',32,0
nl	db	10,0



prova:
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
		
		; esempio: stampa input->e
		mov EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
        ; [EAX]      input->seq; 			// sequenza
		; [EAX + 4]  input->N; 			    // numero elementi sequenza
		; [EAX + 8]  input->sd;			    // seed
		; [EAX + 12] input->to;			    // temperatura
		; [EAX + 16] input->alpha;		    // tasso raffredamento
		; [EAX + 20] input->k; 		        // costante
		; [EAX + 24] input->hydrophobicity;	// hydrophobicity
		; [EAX + 28] input->volume;		    // volume
		; [EAX + 32] input->charge;		    // charge
		; [EAX + 36] input->phi;		    // vettore angoli phi
		; [EAX + 40] input->psi;		    // vettore angoli psi
		; [EAX + 44] input->e;			    // energy
		; [EAX + 48] input->dispaly;
		; [EAX + 52] input->silent;
		MOVSS XMM0, [EAX+44]
		MOVSS [e], XMM0 
		prints msg            
		printss e   
		prints nl
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
        MOV EBX, [EBP+S]
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