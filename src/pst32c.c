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
#define REAL_MATRIX	type**
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
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
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

REAL_MATRIX alloc_real_matrix(int n, int m){
	type** mat; 
	int i, k;

	mat= (type**) _mm_malloc(n*sizeof(type*),16);
	
	for(i=0; i<n; i++){
		mat[i]=(type*) _mm_malloc(m*sizeof(type),16);
		
		if(mat[i]==NULL){
			for(k=0; k<i; k++){ free(mat[k]); }
		}
	}
	return mat;
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

void sub2(VECTOR v1, VECTOR v2, VECTOR res){	
    	for(int i=0;i< 4; i++){
          	res[i]=v1[i]-v2[i];
	}
}

void div_scalareC(VECTOR v, type s){
	if(s==0) s=1;
    for(int i=0;i<3; i++){
        v[i]=v[i]/s;
	}
}

void sum2(VECTOR v1, VECTOR v2, VECTOR res){
	for(int i=0;i< 4; i++){
        res[i]=v1[i]+v2[i];
	}
}

extern void euclidean_dist_sse(VECTOR v1, VECTOR v2, type* res);

void prodotto_vettore_scalare(VECTOR v, type s){
	for(int i=0;i<4;i++){
		v[i]=v[i] * s;
	}
}

void prodotto_scalare2(VECTOR v1, VECTOR v2, type* res){
	type prodotto = 0;
	for(int i=0; i<3;i++){
		prodotto += v1[i]*v2[i];
	}
	*res = prodotto;
}


void prodotto_vettore_matrice2(VECTOR v, REAL_MATRIX m, VECTOR res){
    for (int i = 0; i < 3; i++) { 
        for (int j = 0; j < 3; j++) { 
            res[i] += v[j] * m[j][i];
        }
    }
}

void mul2(VECTOR v1, VECTOR v2){
	for(int i=0;i<4;i++){
		v1[i] = v1[i] * v2[i];
	}
}

type approx_cos(type x) {
    type x2 = x * x;
    return 1 - (x2 / 2.0f) + (x2 * x2 / 24.0f) - (x2 * x2 * x2 / 720.0f);
}

type approx_sin(type x) {
    type x2 = x * x;
    return x - (x * x2 / 6.0f) + (x * x2 * x2 / 120.0f) - (x * x2 * x2 * x2 / 5040.0f);
}



void mul_matrix(type* vector, type** matrix, type* result, int n){
    for (int i = 0; i < n; i++) {
        result[i] = 0.0f;
    }
    
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n; j++) { 
            result[i] += vector[j] * matrix[j][i];
        }
    }
}


void rotation(VECTOR axis, type theta, REAL_MATRIX rotation_matrix){
	type res;
	type a, b, c, d;
	prodotto_scalare2(axis, axis, &res);

	div_scalareC(axis, res);
	prodotto_vettore_scalare(axis, approx_sin(theta / 2.0f));
	prodotto_vettore_scalare(axis, -1);

	a=approx_cos(theta / 2.0f);
	b=axis[0];
    c=axis[1];
	d=axis[2];

	rotation_matrix[0][0] = a * a + b * b - c * c - d * d;
    rotation_matrix[0][1] = 2 * (b * c + a * d);
   	rotation_matrix[0][2] = 2 * (b * d - a * c);

   	rotation_matrix[1][0] = 2 * (b * c - a * d);
	rotation_matrix[1][1] = a * a + c * c - b * b - d * d;
   	rotation_matrix[1][2] = 2 * (c * d + a * b);

  	rotation_matrix[2][0] = 2 * (b * d + a * c);
	rotation_matrix[2][1] = 2 * (c * d - a * b);
	rotation_matrix[2][2] = a * a + d * d - b * b - c * c;
}


MATRIX backbone(char* s, int n, VECTOR phi, VECTOR psi){

	type r_ca_n = 1.46;
	type r_ca_c = 1.52;
	type r_c_n = 1.33;
	
	type theta_ca_c_n = 2.028;
	type theta_c_n_ca = 2.124;
	type theta_n_ca_c = 1.940;
	
	MATRIX coords = alloc_matrix( n*3, 3);
	
	int i, j;
	
	coords[0]=0;
	coords[1]=0;
	coords[2]=0;

	coords[3] = r_ca_n;
	coords[4] = 0; 
	coords[5] = 0;

	VECTOR coords_c = alloc_matrix(1, 4);
	VECTOR coords_c_alpha = alloc_matrix(1, 4);
	VECTOR coords_n = alloc_matrix(1,4);
	VECTOR v1=alloc_matrix(1,4);
	VECTOR v2 = alloc_matrix(1,4);
	VECTOR v3 = alloc_matrix(1,4);

	REAL_MATRIX rotation_matrix = alloc_real_matrix(3,3);
	VECTOR v = alloc_matrix(1,3);
	VECTOR v_ = alloc_matrix(1,3);
	VECTOR newv = alloc_matrix(1,4);

	type norm;

	for(i=0; i<n; i++){
		
		int idx=i*3*3; 
		
		if (i>0) {

			coords_c_alpha[0] = coords[(i-1)*3*3+3];
			coords_c_alpha[1] = coords[(i-1)*3*3+4];
			coords_c_alpha[2] = coords[(i-1)*3*3+5];
			coords_c_alpha[3] = 0;

			coords_c[0] = coords[(i-1)*3*3+6]; 
			coords_c[1] = coords[(i-1)*3*3+7];
			coords_c[2] = coords[(i-1)*3*3+8];
			coords_c[3] = 0;

			v1[0] = coords_c[0]-coords_c_alpha[0];
			v1[1] = coords_c[1]-coords_c_alpha[1];
			v1[2] = coords_c[2]-coords_c_alpha[2];
			v1[3] = coords_c[3]-coords_c_alpha[3];
			
			norm = sqrtf(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] + v1[3]*v1[3]);
			v1[0] = v1[0]/norm;
			v1[1] = v1[1]/norm;
			v1[2] = v1[2]/norm;
			v1[3] = v1[3]/norm;
			
			rotation(v1, theta_c_n_ca, rotation_matrix);
			
			v[0]=0;
			v[1]=r_c_n;
			v[2]=0;
			
			newv[0]=0;
			newv[1]=0;
			newv[2]=0;
			newv[3]=0;

			prodotto_vettore_matrice2(v, rotation_matrix, newv);
			
			coords_n[0] = coords_c[0]+newv[0];
			coords_n[1] = coords_c[1]+newv[1];
			coords_n[2] = coords_c[2]+newv[2];
			coords_n[3] = coords_c[3]+newv[3];
			
			coords[idx] = coords_n[0];
			coords[idx+1] = coords_n[1];
			coords[idx+2] = coords_n[2];

			v2[0] = coords_n[0]-coords_c[0];
			v2[1] = coords_n[1]-coords_c[1];
			v2[2] = coords_n[2]-coords_c[2];
			v2[3] = coords_n[3]-coords_c[3];
			
			norm = sqrtf(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2] + v2[3]*v2[3]);
			v2[0] = v2[0]/norm;
			v2[1] = v2[1]/norm;
			v2[2] = v2[2]/norm;
			v2[3] = v2[3]/norm;
			
			rotation(v2, phi[i], rotation_matrix);
			
			v_[0]=0;
			v_[1]=r_ca_n;
			v_[2]=0;

			newv[0]=0;
			newv[1]=0;
			newv[2]=0;
			newv[3]=0;

			prodotto_vettore_matrice2(v_, rotation_matrix, newv);
			
			coords_c_alpha[0] = coords_n[0]+newv[0];
			coords_c_alpha[1] = coords_n[1]+newv[1];
			coords_c_alpha[2] = coords_n[2]+newv[2];
			coords_c_alpha[3] = coords_n[3]+newv[3];

			coords[idx+3] = coords_c_alpha[0];
			coords[idx+4] = coords_c_alpha[1];
			coords[idx+5] = coords_c_alpha[2];
		}
	
		coords_n[0] = coords[idx];
		coords_n[1] = coords[idx+1];
		coords_n[2] = coords[idx+2];
		coords_n[3] = 0;

		coords_c_alpha[0] = coords[idx+3];
		coords_c_alpha[1] = coords[idx+4];
		coords_c_alpha[2] = coords[idx+5];
		coords_c_alpha[3] = 0; 
	
		v3[0] = coords_c_alpha[0] - coords_n[0];
		v3[1] = coords_c_alpha[1] - coords_n[1];
		v3[2] = coords_c_alpha[2] - coords_n[2];
		v3[3] = coords_c_alpha[3] - coords_n[3];
	
		norm = sqrtf(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] + v3[3]*v3[3]);
		v3[0] = v3[0]/norm;
		v3[1] = v3[1]/norm;
		v3[2] = v3[2]/norm;
		v3[3] = v3[3]/norm;

		rotation(v3, psi[i], rotation_matrix);
		
		newv[0]=0;
		newv[1]=0;
		newv[2]=0;
		newv[3]=0;
	
		v[0] = 0; 
		v[1] = r_ca_c; 
		v[2] = 0;
		
		prodotto_vettore_matrice2(v, rotation_matrix, newv);
	
		coords_c[0] = coords_c_alpha[0] +newv[0];
		coords_c[1] = coords_c_alpha[1] +newv[1];
		coords_c[2] = coords_c_alpha[2] +newv[2];
		coords_c[3] = coords_c_alpha[3] +newv[3];
	
		coords[idx+6] = coords_c[0];
		coords[idx+7] = coords_c[1];
		coords[idx+8] = coords_c[2];
	}

	return coords;
}

type rama_energy(VECTOR phi, VECTOR psi){
	int n = 256;
	type alpha_phi = -57.8;
	type alpha_psi = -47.0;
	type beta_phi = -119.0;
	type beta_psi = 113.0;
	type energy = 0;

	int i;
    
		type min, alpha_dist, beta_dist;
        for(i=0; i<n; i+=4){
			alpha_dist = sqrtf((phi[i]-alpha_phi)*(phi[i]-alpha_phi)+((psi[i]-alpha_psi)*(psi[i]-alpha_psi)));
			beta_dist = sqrtf((phi[i]-beta_phi)*(phi[i]-beta_phi)+(psi[i]-beta_psi)*(psi[i]-beta_psi));
			min = fmin(alpha_dist, beta_dist);
			energy = energy + 0.5 * min;

			alpha_dist = sqrtf((phi[i+1]-alpha_phi)*(phi[i+1]-alpha_phi)+((psi[i+1]-alpha_psi)*(psi[i+1]-alpha_psi)));
			beta_dist = sqrtf((phi[i+1]-beta_phi)*(phi[i+1]-beta_phi)+(psi[i+1]-beta_psi)*(psi[i+1]-beta_psi));
			min = fmin(alpha_dist, beta_dist);
			energy = energy + 0.5 * min;

			alpha_dist = sqrtf((phi[i+2]-alpha_phi)*(phi[i+2]-alpha_phi)+((psi[i+2]-alpha_psi)*(psi[i+2]-alpha_psi)));
			beta_dist = sqrtf((phi[i+2]-beta_phi)*(phi[i+2]-beta_phi)+(psi[i+2]-beta_psi)*(psi[i+2]-beta_psi));
			min = fmin(alpha_dist, beta_dist);
			energy = energy + 0.5 * min;

			alpha_dist = sqrtf((phi[i+3]-alpha_phi)*(phi[i+3]-alpha_phi)+((psi[i+3]-alpha_psi)*(psi[i+3]-alpha_psi)));
			beta_dist = sqrtf((phi[i+3]-beta_phi)*(phi[i+3]-beta_phi)+(psi[i+3]-beta_psi)*(psi[i+3]-beta_psi));
			min = fmin(alpha_dist, beta_dist);
			energy = energy + 0.5 * min;
		
		}
	return energy;
}


type euclidean_dist(VECTOR v1, VECTOR v2, type* res){
	*res = sqrtf((v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
}



type hydrophobic_energy(char* s, int n, MATRIX coords){
	type energy = 0;
	int i,j,k;
	VECTOR coords_c_alpha_i = alloc_matrix(1,4);
	VECTOR coords_c_alpha_j = alloc_matrix(1,4);

	for(i=0; i<n ; i++){
		
		int idx_i = i*3*3+3;
	
		coords_c_alpha_i[0] = coords[idx_i];
		coords_c_alpha_i[1] = coords[idx_i+1];
		coords_c_alpha_i[2] = coords[idx_i+2];
		coords_c_alpha_i[3] = 0;

		for(j=i+1; j<n ; j++){

			int idx_j=j*3*3+3; 

			coords_c_alpha_j[0] = coords[idx_j];
			coords_c_alpha_j[1] = coords[idx_j+1];
			coords_c_alpha_j[2] = coords[idx_j+2];
			coords_c_alpha_j[3] = 0;

			type dist = 0;

			euclidean_dist_sse(coords_c_alpha_i, coords_c_alpha_j, &dist);

			int pos_i = s[i]-65;
            int pos_j = s[j]-65;

			if(dist < 10.0)
				energy = energy + (hydrophobicity[pos_i]*hydrophobicity[pos_j])/dist;
		}
	}
	return energy;
}


type electrostatic_energy(char* s, int n, MATRIX coords){
        type energy = 0;
        int i,j,k;
		VECTOR coords_c_alpha_i = alloc_matrix(1,4);
		VECTOR coords_c_alpha_j = alloc_matrix(1,4);

        for(i=0; i<n ; i++){

			int idx_i = i*3*3+3;
	
			coords_c_alpha_i[0] = coords[idx_i];
			coords_c_alpha_i[1] = coords[idx_i+1];
			coords_c_alpha_i[2] = coords[idx_i+2];
			coords_c_alpha_i[3] = 0;

			for(j=i+1; j<n ; j++){
			
				int idx_j = j*3*3+3;
					
				coords_c_alpha_j[0] = coords[idx_j];
				coords_c_alpha_j[1] = coords[idx_j+1];
				coords_c_alpha_j[2] = coords[idx_j+2];
				coords_c_alpha_j[3] = 0;
				type dist = 0;

				euclidean_dist_sse(coords_c_alpha_i, coords_c_alpha_j, &dist);
				int pos_i = s[i]-65;
				int pos_j = s[j]-65;

				if(i!=j && dist < 10.0 && charge[pos_i]!=0 && charge[pos_j]!=0)
					energy = energy + (charge[pos_i]*charge[pos_j])/(dist*4.0);
					
			}
        return energy;
	}
}


type packing_energy(char* s, int n, MATRIX coords){
        type energy = 0;
        int i,j,k;
		VECTOR coords_c_alpha_i = alloc_matrix(1,4);
		VECTOR coords_c_alpha_j = alloc_matrix(1,4);

        for(i=0; i<n ; i++){
            type density = 0;
            int pos_i = s[i]-65;

            int idx_i = i*3*3+3; 
                    
            coords_c_alpha_i[0] = coords[idx_i];
			coords_c_alpha_i[1] = coords[idx_i+1];
			coords_c_alpha_i[2] = coords[idx_i+2];
			coords_c_alpha_i[3] = 0;
                
		    for(j=0; j<n ; j++){

                int idx_j=j*3*3+3; 
		
                coords_c_alpha_j[0] = coords[idx_j];
				coords_c_alpha_j[1] = coords[idx_j+1];
				coords_c_alpha_j[2] = coords[idx_j+2];
				coords_c_alpha_j[3] = 0;
                       
		       	type dist = 0;
				euclidean_dist_sse(coords_c_alpha_i, coords_c_alpha_j, &dist);
                int pos_j = s[j]-65;

                if(i!=j && dist < 10.0)
                    density = density + (volume[pos_j])/(dist*dist*dist);
		}
		energy = energy +((volume[pos_i]-density)*(volume[pos_i]-density));
        }
        return energy;
}


type energy(char* s, int n, VECTOR phi, VECTOR psi){
        MATRIX coords = backbone(s, n, phi, psi);

        type rama = rama_energy(phi, psi);
        type hydro = hydrophobic_energy(s,n, coords);
        type elec = electrostatic_energy(s,n, coords);
        type pack = packing_energy(s,n, coords);

		VECTOR v1 = alloc_matrix(1,4);

		v1[0] = rama; v1[1] = hydro; v1[2] = elec; v1[3] = pack;

        type w_rama = 1.0;
        type w_hydro = 0.5;
        type w_elec = 0.2;
        type w_pack = 0.3;

		VECTOR v2 = alloc_matrix(1, 4);

       	v2[0] = w_rama; v2[1] = w_hydro; v2[2] = w_elec; v2[3] = w_pack;
	
        mul2(v1, v2);
		type total=0;
		for(int i=0; i<4 ; i++) 
			total+=v1[i];

		return total;
}


void pst(params* input){
	// --------------------------------------------------------------
	// Codificare qui l'algoritmo di Predizione struttura terziaria
	// --------------------------------------------------------------
	type E = energy(input->seq,input->N, input->phi, input->psi);
	type T = input->to;

	int t=0;
	printf("Energia iniziale:%f\n",E);

	while(T > 0){
		unsigned int seed = input->sd;
		int i = (int) (random()*input->N); 

		type theta_phi = (random()*2*M_PI)-M_PI;
		type theta_psi = (random()*2*M_PI)-M_PI;

		input->phi[i] = input->phi[i] + theta_phi;
		input->psi[i] = input->psi[i] + theta_psi;

		type new_E = energy(input->seq,input->N, input->phi, input->psi);
		type delta_E = new_E - E;

		if(delta_E <= 0){
			E = new_E;
		}else {
			type P = exp(-delta_E / (input->k * T));
			type r = (type) rand() / RAND_MAX; 

			if(r<=P){
				E = new_E;
			}else{
				input->phi[i] = input->phi[i] - theta_phi;
				input->psi[i] = input->psi[i] - theta_psi;
			}
		}
		t = t+1;
		T = input->to - sqrt((input->alpha)*t);
	}
	input->e=E;

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
