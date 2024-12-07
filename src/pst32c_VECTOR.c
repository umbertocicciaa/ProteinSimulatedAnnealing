/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2024/25
* 
* Progetto dell'algoritmo Predizione struttura terziaria proteine 221 231 a
* in linguaggio assembly x86-32 + SSE
* 
* F. Angiulli F. Fassetti S. Nisticò, novembre 2024
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm && ./pst32c $pars
* 
* oppure
* 
* ./runpst32
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <malloc.h>
#include <stdint.h>

#define	type		float
#define	MATRIX		type*
#define REAL_MATRIX type**
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge

typedef struct {
	char* seq;		// sequenza di amminoacidi
	int N;			// lunghezza sequenza
	unsigned int sd; 	// seed per la generazione casuale
	type to;		// temperatura INIZIALE
	type alpha;		// tasso di raffredamento
	type k;		// costante
	VECTOR hydrophobicity; // hydrophobicity
	VECTOR volume;		// volume
	VECTOR charge;		// charge
	VECTOR phi;		// vettore angoli phi
	VECTOR psi;		// vettore angoli psi
	type e;		// energy
	int display;
	int silent;

} params;


/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) {  
    if (size <= 0 || elements <= 0 || size > SIZE_MAX / elements) {
        return NULL; // Dimensioni non valide o overflow
    }
    return _mm_malloc(elements * size, 16); 
}

void free_block(void* p) { 
	_mm_free(p);
}

REAL_MATRIX alloc_real_matrix(int n, int m) {
    REAL_MATRIX mat;
    int i;

    // Alloca la memoria per un array di puntatori (con allineamento a 16 byte)
    mat = (REAL_MATRIX) malloc(n * sizeof(type*));
	
    // Alloca la memoria per ciascuna riga (con allineamento a 16 byte)
    for (i = 0; i < n; i++) {
        mat[i] = (type*) malloc(m * sizeof(type));
    }

    return mat;  // Restituisci la matrice allocata
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}


int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}



/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 
* 	load_seq
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	char* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(char), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero floating-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
* 	gen_rnd_mat
* 	=========
* 
*	Genera in maniera casuale numeri reali tra -pi e pi
*	per riempire una struttura dati di dimensione Nx1
* 
*/
void gen_rnd_mat(VECTOR v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random()*2 * M_PI) - M_PI;
	}
}

// PROCEDURE ASSEMBLY

//extern void prova(params* input);
extern void sum(VECTOR v1, VECTOR v2, VECTOR res);
extern void sub(VECTOR v1, VECTOR v2, VECTOR res);
extern void euclidean_dist_sse(VECTOR v1, VECTOR v2, type* res);

void sottrazione(VECTOR v1, VECTOR v2, VECTOR res)
{	
    for(int i=0;i< 3; i++){
           res[i]=v1[i]-v2[i];
		   
	}
}

void div_scalare(VECTOR v, type* s){
	if(*s==0){
		for(int i=0;i<3; i++){
        	v[i]=v[i]/(1);
		}	
	}else{
		for(int i=0;i<3; i++){
        	v[i]=v[i]/(*s);
		}
	}
	
}

void somma(VECTOR v1, VECTOR v2, VECTOR res){
	for(int i=0;i< 3; i++){
        res[i]=v1[i]+v2[i];
	}
}


void norm(VECTOR v, type* res){
	type norm = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	
	norm = sqrtf(norm);
	res = &norm;

}

void prodotto_vettore_scalare(VECTOR v, type s){
	for(int i=0;i<3;i++){
		v[i]=v[i] * s;
	}
}

void prodotto_scalare(VECTOR v1, VECTOR v2, type* res){
	type prodotto = 0;
	for(int i=0; i<3;i++){
		prodotto += v1[i]*v2[i];
	}
	*res = prodotto;
}

/*void prodotto_vettore_matrice(VECTOR v, REAL_MATRIX m, VECTOR res){
	for (int i = 0; i < 3; i++) {
        res[i] = 0.0f;
    }
    
    // Esegui il prodotto tra il vettore e la matrice
    for (int i = 0; i < 3; i++) { // Itera sulle colonne della matrice
        for (int j = 0; j < 3; j++) { // Itera sulle righe della matrice
            res[i] += v[j] * m[j][i];
			printf("%f += %f*%f\n", res[i], v[j], m[j][i]);
        }
    }
}*/

void prodotto_vettore_matrice(VECTOR v, MATRIX m, VECTOR res){
    // Esegui il prodotto tra il vettore e la matrice
    for (int i = 0; i < 3; i++) { // Itera sulle colonne della matrice
        for (int j = 0; j < 3; j++) { // Itera sulle righe della matrice
            res[i] += v[j] * m[i+j];
        }
    }
}

void mul(VECTOR v1, VECTOR v2){
	for(int i=0;i<3;i++){
		v1[i] = v1[i] * v2[i];
	}
}

// Funzioni per approssimazioni polinomiali di coseno e seno
type approx_cos(type x) {
    type x2 = x * x;
    return 1 - (x2 / 2.0f) + (x2 * x2 / 24.0f) - (x2 * x2 * x2 / 720.0f);
}

type approx_sin(type x) {
    type x2 = x * x;
    return x - (x * x2 / 6.0f) + (x * x2 * x2 / 120.0f) - (x * x2 * x2 * x2 / 5040.0f);
}

void normalize_axis(VECTOR axis) {
    // Calcola la norma del vettore axis (prodotto scalare di axis con se stesso)
    type norm = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    norm = sqrtf(norm); // Calcolo della radice quadrata
    // Normalizza il vettore axis
	if(norm==0) norm=1;
    axis[0] /= norm;
    axis[1] /= norm;
    axis[2] /= norm;
}


void mul_matrix(type* vector, type** matrix, type* result, int n){
	// Inizializza il risultato a zero
    for (int i = 0; i < n; i++) {
        result[i] = 0.0f;
    }
    
    // Esegui il prodotto tra il vettore e la matrice
    for (int i = 0; i < n; i++) { // Itera sulle colonne della matrice
        for (int j = 0; j < n; j++) { // Itera sulle righe della matrice
            result[i] += vector[j] * matrix[j][i];
        }
    }
}


/*void rotation(VECTOR axis, type theta, REAL_MATRIX rotation_matrix) {
    // Normalizzazione dell'asse
    normalize_axis(axis);

    // Calcolo di a, b, c, d
    type a = approx_cos(theta / 2.0f);
    type b = -axis[0] * approx_sin(theta / 2.0f);
    type c = -axis[1] * approx_sin(theta / 2.0f);
    type d = -axis[2] * approx_sin(theta / 2.0f);

    // Calcolo dei termini della matrice di rotazione
    rotation_matrix[0][0] = a * a + b * b - c * c - d * d;
    rotation_matrix[0][1] = 2 * (b * c + a * d);
    rotation_matrix[0][2] = 2 * (b * d - a * c);

    rotation_matrix[1][0] = 2 * (b * c - a * d);
    rotation_matrix[1][1] = a * a + c * c - b * b - d * d;
    rotation_matrix[1][2] = 2 * (c * d + a * b);

    rotation_matrix[2][0] = 2 * (b * d + a * c);
    rotation_matrix[2][1] = 2 * (c * d - a * b);
    rotation_matrix[2][2] = a * a + d * d - b * b - c * c;
}*/

/*void rotation(VECTOR axis, type theta, REAL_MATRIX rotation_matrix){
	type* res;
	
	type a, b, c, d;
	prodotto_scalare(axis, axis, res);

	printf("axis[");

		for (int i = 0; i < 3; i++) {
			printf("%f, ", axis[i]);
		}
		printf("]\n");

	type val = *res;
	div_scalare(axis, val);
	
	printf("axis DIV SCALARE[");

		for (int i = 0; i < 3; i++) {
			printf("%f, ", axis[i]);
		}
		printf("]\n");
	prodotto_vettore_scalare(axis, approx_sin(theta / 2.0f));
	prodotto_vettore_scalare(axis, -1);


	printf("axis VETORE SCALARE[");

		for (int i = 0; i < 3; i++) {
			printf("%f, ", axis[i]);
		}
		printf("]\n");

	a=approx_cos(theta / 2.0f);
	b=axis[0];
    c=axis[1];
	d=axis[2];

	printf("a: %f, b: %f, c: %f, d: %f\n", a, b, c,d);

	// Calcolo dei termini della matrice di rotazione
	rotation_matrix[0][0] = a * a + b * b - c * c - d * d;
    rotation_matrix[0][1] = 2 * (b * c + a * d);
   	rotation_matrix[0][2] = 2 * (b * d - a * c);

   	rotation_matrix[1][0] = 2 * (b * c - a * d);
	rotation_matrix[1][1] = a * a + c * c - b * b - d * d;
   	rotation_matrix[1][2] = 2 * (c * d + a * b);

  	rotation_matrix[2][0] = 2 * (b * d + a * c);
	rotation_matrix[2][1] = 2 * (c * d - a * b);
	rotation_matrix[2][2] = a * a + d * d - b * b - c * c;

	for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.2f ", rotation_matrix[i][j]);
        }
        printf("\n");
	}
}*/

void rotation(VECTOR axis, type theta, MATRIX rotation_matrix){
	type res=0;
	type a, b, c, d;

	prodotto_scalare(axis, axis, &res);

	div_scalare(axis, &res);

	prodotto_vettore_scalare(axis, approx_sin(theta / 2.0f));
	
	prodotto_vettore_scalare(axis, -1);

	a=approx_cos(theta / 2.0f);
	b=axis[0];
    c=axis[1];
	d=axis[2];


	// Calcolo dei termini della matrice di rotazione
	rotation_matrix[0] = a * a + b * b - c * c - d * d;
    rotation_matrix[1] = 2 * (b * c + a * d);
   	rotation_matrix[2] = 2 * (b * d - a * c);

   	rotation_matrix[3] = 2 * (b * c - a * d);
	rotation_matrix[4] = a * a + c * c - b * b - d * d;
   	rotation_matrix[5] = 2 * (c * d + a * b);

  	rotation_matrix[6] = 2 * (b * d + a * c);
	rotation_matrix[7] = 2 * (c * d - a * b);
	rotation_matrix[8] = a * a + d * d - b * b - c * c;
	

}


MATRIX backbone(int n, VECTOR phi, VECTOR psi){
	
	//distanze atomi in Angstrom
	type r_ca_n = 1.46;
	type r_ca_c = 1.52;
	type r_c_n = 1.33;
	
	//angoli atomi in radianti
	type theta_ca_c_n = 2.028;
	type theta_c_n_ca = 2.124;
	type theta_n_ca_c = 1.940;
	
	//allocazione matrice coords
	//MATRIX coords = alloc_matrix( n*3, 3);
	MATRIX coords = alloc_matrix( 1, n*3*3);

	int i, j;
	
	//inizializzazione coordinate N primo amminoacido
	for(i=0; i<3; i++) 
		coords[i]=0;
	
	//iniziaizzazione coordinate C_alpha primo amminoacido
	i=3;
	coords[i] = r_ca_n;
	coords[i+1] = 0; 
	coords[i+2] = 0; //allocato fino a coords[5]
	
	for(i=0; i<n; i++){
		
		int idx=i*3*3; //calcolo indice base amminoacido
		
		
		if (i > 0) {
			//posiziona N
			VECTOR coords_c = alloc_matrix(1, 4);
			VECTOR coords_c_alpha = alloc_matrix(1, 4);
			for(j=0; j<3; j++){
				coords_c_alpha[j] = coords[(i-1)*3*3+j+3];
				coords_c[j] = coords[(i-1)*3*3+j+6];
			}

			j=3;
			coords_c_alpha[j]=0;
			coords_c[j]=0;
		
	
			VECTOR v1 = alloc_matrix(1, 4);

			// C-C_alpha
			sub(coords_c, coords_c_alpha, v1);

			normalize_axis(v1);

			//costruzione matrice rotation
			//REAL_MATRIX rotation_matrix = alloc_real_matrix(3,3);
			MATRIX rotation_matrix = alloc_matrix(3,3);
			rotation(v1, theta_c_n_ca, rotation_matrix);

			VECTOR v = alloc_matrix(1,3);
			v[0]=0.0f;
			v[1]=r_c_n;
			v[2]=0.0f;
			

			VECTOR newv = alloc_matrix(1,4);
			prodotto_vettore_matrice(v, rotation_matrix, newv);
		
			VECTOR coords_n = alloc_matrix(1,4);

			//N = C + newv
			sum(coords_c, newv, coords_n);
			

			//inserisco in coords le coordinate di N
			for(j=0; j<3; j++)	coords[idx+j] = coords_n[j];
			
			//posiziona C_alpha
			VECTOR v2 = alloc_matrix(1,4);

			//N - C
			sub(coords_n, coords_c, v2);

			normalize_axis(v2);
			rotation(v2, phi[i], rotation_matrix);


			VECTOR v_ = alloc_matrix(1,3);
			v_[0]=0.0f;
			v_[1]=r_ca_n;
			v_[2]=0.0f;
			
			prodotto_vettore_matrice(v_, rotation_matrix, newv);

			//C_alpha = N + newv
			sum(coords_n, newv, coords_c_alpha);
			
			//inserisco in coords le coordinate di C_alpha
			for(j=0; j<3; j++){
				coords[idx+3+j] = coords_c_alpha[j];
			}	
			
		}
		//Posiziona C
	
		VECTOR coords_n = alloc_matrix(1,4);
		VECTOR coords_c_alpha = alloc_matrix(1,4);
		
		VECTOR v3 = alloc_matrix(1,4);
		
		for(int j=0; j<3; j++){
			coords_n[j] = coords[idx+j];
			coords_c_alpha[j] = coords[idx+3+j];
		}

		j=3;
		coords_c_alpha[j]=0;
		coords_n[j]=0;


		//C_alpha - N
		sub(coords_n, coords_c_alpha, v3); //chiamata Assembly

		normalize_axis(v3);

		MATRIX rotation_matrix = alloc_matrix(3,3);
		
		rotation(v3, psi[i], rotation_matrix);		

		VECTOR newv = alloc_matrix(1,3);
		VECTOR v = alloc_matrix(1,3);
		
		v[0]=0;
		v[1]=r_ca_c;
		v[2]=0;


		prodotto_vettore_matrice(v, rotation_matrix, newv);

		//C = C_alpha + newv
		VECTOR coords_c = alloc_matrix(1,4);
		sum(coords_c_alpha, newv, coords_c);
		

		//inserisco in coords le coordinate di C
		for(int j=0; j<3; j++)    coords[idx+6+j] = coords_c[j];
	}
	return coords;
}

/*REAL_MATRIX backbone(char* s, VECTOR phi, VECTOR psi){
	int n =strlen(s);	

	//distanze standard nel backbone
	type r_CA_N = 1.46;
	type r_CA_C = 1.52;
	type r_C_N = 1.33;

	//angoli standard in radianti nel backbone
	type theta_CA_C_N = 2.028;
	type theta_C_N_CA = 2.124;
	type theta_N_CA_C = 1.940;

	//crea la matrice coords 
	REAL_MATRIX coords = alloc_real_matrix(n,9);
	//MATRIX coords = alloc_matrix(n*3, 3);
	coords[0][0] = 0;//x di N nel primo amminoacido
	coords[0][1] = 0;//y di N nel primo amminoacido
	coords[0][2] = 0;//z di N nel primo amminoacido 
	
	coords[0][3] = r_CA_N; //x di Ca nel primo amminoacido
	coords[0][4] = 0; //y di Ca nel primo amminoacido
	coords[0][5] = 0; //z di Ca nel primo amminoacido
	
	int i;
	for(i=0; i<n; i++){
		if(i>0){
			VECTOR v1 = alloc_matrix(1,3);

			int j;
			for(j=0; j<3; j++){
				v1[j] = coords[i-1][j+6] - coords[i-1][j+3]; //differenza C-Ca
			}
			normalize_axis(v1);
			REAL_MATRIX rot_v1 = alloc_real_matrix(3,3);
			rotation(v1, theta_C_N_CA, rot_v1);
			
			type v_v1[3] = {0.0, r_C_N, 0.0};

			VECTOR newv_v1 = alloc_matrix(1,3);
			mul_matrix(v_v1, rot_v1, newv_v1, 3);

			for(j=0; j<3; j++){
				coords[i][j] = coords[i-1][j+6] + newv_v1[j]; //posiziona N usando newv e C dell'amminoacido precedente
			}

			//POSIZIONA Ca USANDO PHI
			VECTOR v2 = alloc_matrix(1,3);

			for(j=0; j<3;j++){
				v2[j] = coords[i][j] - coords[i-1][j+6];
			}
			normalize_axis(v2);
			REAL_MATRIX rot_v2 = alloc_real_matrix(3,3);
			rotation(v2, phi[i], rot_v2);
			type v_v2[3] = {0, r_CA_N, 0};
			VECTOR newv_v2 = alloc_matrix(1,3);
			mul_matrix(v_v2, rot_v2, newv_v2, 3);

			for(j=0; j<3;j++){
				coords[i][j+3] = coords[i][j] + newv_v2[j]; //posiziona Ca usando N dell'amminoacido corrente
			}
		}

		//posiziona C usando psi
		VECTOR v3 = alloc_matrix(1,3);

		int j;
		for(j=0; j<3; j++){
			v3[j] = coords[i][j+3] - coords[i][j];
			// coords[0][6] = coords[0][3] - coords[0][0];
			// coords[0][7] = coords[0][4] - coords[0][1];
			// coords[0][8] = coords[0][5] - coords[0][2];
		}

		normalize_axis(v3);
		
		REAL_MATRIX rot = alloc_real_matrix(3,3);
		rotation(v3, psi[i],rot);
		
		type v[3] = {0, r_CA_C, 0};
		VECTOR newv = alloc_matrix(1,3);
		mul_matrix(v, rot, newv, 3);

		for(j=0; j<3; j++){
			coords[i][j+6] = newv[j] + coords[i][j+3]; //aggiorna C
		}
	}
	return coords;
}*/


type rama_energy(VECTOR phi, VECTOR psi){
	int n = 26;
	type alpha_phi = -57.8;
	type alpha_psi = -47.0;
	type beta_phi = -119.0;
	type beta_psi = 113.0;
	type energy = 0;

	int i;
	for(i=0; i<n; i++){
		type alpha_dist = sqrt(pow(phi[i]-alpha_phi,2)+pow(psi[i]-alpha_psi,2));
		type beta_dist = sqrt(pow(phi[i]-beta_phi,2)+pow(psi[i]-beta_psi,2));
		type min;
		if(alpha_dist<beta_dist) min = alpha_dist;
		else min = beta_dist;
		energy = energy + 0.5 * min;
	}
	return energy;
}

/*type euclidean_distance(VECTOR v1, VECTOR v2){
	type res = sqrtf(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2) + pow(v2[2]-v1[2],2));
	return res;
}*/


/*type hydrophobic_energy(char* s, REAL_MATRIX coords){
	int n = strlen(s);
	type energy = 0;

	int i,j;
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++){
			//considera la distanza euclidea tra gli Ca degli amminoacidi in pos i e j
			type dist = euclidean_dist(coords, i, j);
			int pos_i = s[i]-65;
			int pos_j = s[j]-65;

			if(dist<10.0){
				energy = energy + (hydrophobicity[pos_i]*hydrophobicity[pos_j])/dist;
			}
		}
	}
	return energy;

}*/

type hydrophobic_energy(char* s, int n, MATRIX coords){
	
	type energy = 0;
	int i,j,k;

	//estrapolo coordinate C_alpha
	for(i=0; i<n ; i++){
		
		int idx_i = i*3*3+3; //indirizzo base generico C_alpha_i
		
		VECTOR coords_c_alpha_i = alloc_matrix(1,4);
		for(k=0; k<3 ; k++)	coords_c_alpha_i[k] = coords[idx_i+k];
		k=3;
		coords_c_alpha_i[k]=0;
		
		for(j=i+1; j<n ; j++){

			int idx_j=j*3*3+3; //indirizzo base generico C_alpha_j

			VECTOR coords_c_alpha_j = alloc_matrix(1,4);
			for(k=0; k<3 ; k++) coords_c_alpha_j[k] = coords[idx_j+k];
			k=3;
			coords_c_alpha_j[k]=0;
			
			
			type dist=0;
			euclidean_dist_sse(coords_c_alpha_i, coords_c_alpha_j,&dist);
			int pos_i = s[i]-65;
            int pos_j = s[j]-65;

			if(dist < 10.0){
				if(dist==0) dist=1;
				energy = energy + (hydrophobicity[pos_i]*hydrophobicity[pos_j])/(dist);
			}
		}
	}
	return energy;
}

/*type electrostatic_energy(char* s, REAL_MATRIX coords){
	int n = strlen(s);
	type energy = 0;

	int i,j;
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++){
			//considera la distanza euclidea tra gli Ca degli amminoacidi in pos i e j
			type dist = euclidean_dist(coords, i, j);
			int pos_i = s[i]-65;
			int pos_j = s[j]-65;

			if(i!=j && dist<10.0 && charge[pos_i]!=0 && charge[pos_j]!=0){
				energy = energy + (charge[pos_i]*charge[pos_j])/(dist*4.0);
			}
		}
	}
	return energy;
}*/

type electrostatic_energy(char* s, int n, MATRIX coords){
	
    type energy = 0;
    int i,j,k;
    //estrapolo coordinate C_alpha
    for(i=0; i<n ; i++){
		int idx_i = i*3*3+3; //indirizzo base generico C_alpha_i
                
		VECTOR coords_c_alpha_i = alloc_matrix(1,4);
        for(k=0; k<3 ; k++)     
			coords_c_alpha_i[k] = coords[idx_i+k];
        k=3;
		coords_c_alpha_i[k]=0;

		for(j=i+1; j<n ; j++){
			int idx_j = j*3*3+3;
            VECTOR coords_c_alpha_j = alloc_matrix(1,4);
            for(k=0; k<3 ; k++) coords_c_alpha_j[k] = coords[idx_j+k];
			k=3;
			coords_c_alpha_j[k]=0;
                        
			type dist=0;
			euclidean_dist_sse(coords_c_alpha_i, coords_c_alpha_j, &dist);
            int pos_i = s[i]-65;
            int pos_j = s[j]-65;

            if(i!=j && dist < 10.0 && charge[pos_i]!=0 && charge[pos_j]!=0){
                if(dist==0) dist=1;
				energy = energy + (charge[pos_i]*charge[pos_j])/((dist)*4.0);
			
			}
        }
    }
    return energy;

}

/*type packing_energy(char* s, REAL_MATRIX coords){
	int n = strlen(s);
	type energy = 0;

	int i,j;
	for(i=0; i<n; i++){
		type density = 0;
		int pos_i = s[i]-65;

		for(j=0; j<n; j++){
			//considera la distanza euclidea tra gli Ca degli amminoacidi in pos i e j
			type dist = euclidean_dist(coords, i, j);
			int pos_j = s[j]-65;

			if(i!=j && dist<10.0){
				density = density + (volume[pos_j])/pow(dist, 3);
			}
		}
		energy = energy + pow((volume[pos_i]-density),2);
	}
	return energy;
}*/

type packing_energy(char* s, int n, MATRIX coords){
	
    type energy = 0;
    int i,j,k;
    //estrapolo coordinate C_alpha_i e C_alpha_j
    for(i=0; i<n ; i++){
		type density = 0;
		int pos_i = s[i]-65;

        int idx_i = i*3*3+3; //indirizzo base generico C_alpha_i
                
		VECTOR coords_c_alpha_i = alloc_matrix(1,4);
        for(k=0; k<3 ; k++)     
			coords_c_alpha_i[k] = coords[idx_i+k];
		k=3;
		coords_c_alpha_i[k]=0;
                
		for(j=i+1; j<n ; j++){
            int idx_j=j*3*3+3; //indirizzo base generico C_alpha_j
		
	     	VECTOR coords_c_alpha_j = alloc_matrix(1,4);
            for(k=0; k<3 ; k++) 
				coords_c_alpha_j[k] = coords[idx_j+k];   
			k=3;
			coords_c_alpha_j[k]=0;
			      

		    type dist=0;
			euclidean_dist_sse(coords_c_alpha_i, coords_c_alpha_j, &dist);
			
            int pos_j = s[j]-65;

            if(i!=j && dist < 10.0){
                if(dist==0) dist=1;
				density = density + (volume[pos_j])/pow(dist, 3);
			}
		}
		energy = energy + pow((volume[pos_i]-density),2);
    }
    return energy;
}

/*type energy(char* s, VECTOR phi, VECTOR psi){

	REAL_MATRIX coords = backbone(s, phi, psi);

	//calcolo delle componenti energetiche
	type rama = rama_energy(phi, psi);
	type hydro = hydrophobic_energy(s, coords);
	type elec = electrostatic_energy(s, coords);
	type pack = packing_energy(s, coords);

	//pesi per i diversi contributi
	type w_rama = 1.0;
	type w_hydro = 0.5;
	type w_elec = 0.2;
	type w_pack = 0.3;

	//energia totale
	type total = w_rama*rama + w_hydro*hydro + w_elec*elec + w_pack*pack;
	return total;
}*/

type energy(char* s, int n, VECTOR phi, VECTOR psi){

        MATRIX coords = backbone(n, phi, psi);

        //calcolo delle componenti energetiche
        type rama = rama_energy(phi, psi);
        type hydro = hydrophobic_energy(s, n, coords);
        type elec = electrostatic_energy(s, n, coords);
        type pack = packing_energy(s, n, coords);

		printf("rama: %f, hydro: %f, elec: %f, pack:%f\n",rama,hydro, elec, pack);

		VECTOR v1 = alloc_matrix(1,4);
		//popolo vettore componenti energetiche
		v1[0] = rama; 
		v1[1] = hydro; 
		v1[2] = elec; 
		v1[3] = pack;

        //pesi per i diversi contributi
        type w_rama = 1.0;
        type w_hydro = 0.5;
        type w_elec = 0.2;
        type w_pack = 0.3;

		VECTOR v2 = alloc_matrix(1, 4);
		//popolo vettore dei pesi
       	v2[0] = w_rama; v2[1] = w_hydro; v2[2] = w_elec; v2[3] = w_pack;
	
        //energia totale = w_rama*rama + w_hydro*hydro + w_elec*elec + w_pack*pack; 
        mul(v1, v2);
		type total=0;
		for(int i=0; i<4 ; i++) total+=v1[i];

	return total;
}


void pst(params* input){
	// --------------------------------------------------------------
	// Codificare qui l'algoritmo di Predizione struttura terziaria
	// --------------------------------------------------------------
	

	
	type T = input->to;

	//calcola l'energia
	type E = energy(input->seq, input->N, input->phi, input->psi);
	type t=0;

	while(T > 0){
		//genera un vicino della soluzione corrente
		//int i = rand()%(input->N); //prova anche altra formula col seed!
		unsigned int seed = input->sd;
		int i = rand_r(&seed)%(input->N);
		
		//calcola variazioni casuali
		type theta_phi = (random()*2*M_PI)-M_PI;
		type theta_psi = (random()*2*M_PI)-M_PI;

		//aggiorna i valori degli angoli
		input->phi[i] = input->phi[i] + theta_phi;
		input->psi[i] = input->psi[i] + theta_psi;

		//calcola la variazione dell'energia
		type new_E = energy(input->seq, input->N, input->phi, input->psi);
		type delta_E = new_E - E;

		if(delta_E <= 0){
			//la configurazione è migliorata, la sostituisco a quella vecchia
			E = new_E;
		}else {
			//calcola la probabilità di accettazione
			type P = exp(-delta_E / (input->k * T));
			type r = (type) rand() / RAND_MAX;  

			if(r<=P){
				//accetta la nuova configurazione
				E = new_E;
			}else{
				//rifiuta la configurazione ripristina i valori per psi e phi
				input->phi[i] = input->phi[i] - theta_phi;
				input->psi[i] = input->psi[i] - theta_psi;
			}
		}
		//aggiorna la temperatura
		t = t+1;
		T = input->to - sqrt(input->alpha*t);
	}

}



int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	float time;
	int d;
	
	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->seq = NULL;	
	input->N = -1;			
	input->to = -1;
	input->alpha = -1;
	input->k = -1;		
	input->sd = -1;		
	input->phi = NULL;		
	input->psi = NULL;
	input->silent = 0;
	input->display = 0;
	input->e = -1;
	input->hydrophobicity = hydrophobicity;
	input->volume = volume;
	input->charge = charge;


	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//
	if(argc <= 1){
		printf("%s -seq <SEQ> -to <to> -alpha <alpha> -k <k> -sd <sd> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tSEQ: il nome del file ds2 contenente la sequenza amminoacidica\n");
		printf("\tto: parametro di temperatura\n");
		printf("\talpha: tasso di raffredamento\n");
		printf("\tk: costante\n");
		printf("\tsd: seed per la generazione casuale\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-seq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			seqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-to") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing to value!\n");
				exit(1);
			}
			input->to = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-alpha") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing alpha value!\n");
				exit(1);
			}
			input->alpha = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sd") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing seed value!\n");
				exit(1);
			}
			input->sd = atoi(argv[par]);
			par++;
		}else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(seqfilename == NULL || strlen(seqfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	input->seq = load_seq(seqfilename, &input->N, &d);

	
	if(d != 1){
		printf("Invalid size of sequence file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->to <= 0){
		printf("Invalid value of to parameter!\n");
		exit(1);
	}

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	if(input->alpha <= 0){
		printf("Invalid value of alpha parameter!\n");
		exit(1);
	}

	input->phi = alloc_matrix(input->N, 1);
	input->psi = alloc_matrix(input->N, 1);
	// Impostazione seed 
	srand(input->sd);
	// Inizializzazione dei valori
	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", seqfilename);
		printf("Sequence lenght: %d\n", input->N);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Predizione struttura terziaria
	//
	t = clock();
	pst(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("PST time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname_phi, "out32_%d_%d_phi.ds2", input->N, input->sd);
	save_out(fname_phi, input->phi, input->N);
	sprintf(fname_psi, "out32_%d_%d_psi.ds2", input->N, input->sd);
	save_out(fname_psi, input->psi, input->N);
	if(input->display){
		if(input->phi == NULL || input->psi == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("energy: %f, phi: [", input->e);
			for(i=0; i<input->N; i++){
				printf("%f,", input->phi[i]);
				
			}
			printf("]\n");
			printf("psi: [");
			for(i=0; i<input->N; i++){
				printf("%f,", input->psi[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	
	free(input);

	return 0;
}
