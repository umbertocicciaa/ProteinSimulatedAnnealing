#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <malloc.h>
#include <stdint.h>

#define	type		double
#define	MATRIX		type*
#define REAL_MATRIX	type**
#define	VECTOR		type*



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


int main(){
    int rows;
	int cols;
	//char f[256] = "../src64/phi_double.ds2";
	char f[256] = "./psi_double.ds2";
	MATRIX phi = load_data(f, &rows, &cols);
	printf("rows: %i, cols: %i\n", rows, cols);
	int i;
	printf("psi: [");
	for(i=0; i<rows; i++){
		printf("%f,", phi[i]);
	}
	printf("]\n");
    
}