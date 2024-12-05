#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <malloc.h>
#include <stdint.h>
 
#define type        float
#define MATRIX      type*
#define REAL_MATRIX type**
#define VECTOR      type*
 
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
 
extern void sum(VECTOR v1, VECTOR v2, VECTOR res);
extern void euclidean_dist_sse(VECTOR v1, VECTOR v2, type* res);
 
void somma(VECTOR v1, VECTOR v2, VECTOR res){
    int i;
    for ( i = 0; i < 4; i++){
        res[i]=v1[i]+v2[i];
    }
}
 
type euclidean_distance(VECTOR v1, VECTOR v2){
    type res = sqrtf(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2) + pow(v2[2]-v1[2],2));
    return res;
}
 
int main(){
    VECTOR v1 = alloc_matrix(1,3);
    VECTOR v2 = alloc_matrix(1,3);
    VECTOR v3 = alloc_matrix(1,3);
 
    /*type v1[3]={0.1,2.3,6.7};
    type v2[3]={0.1,2.3,6.7};
    type v3[3]={0.0,0.0,0.0};*/
 
    v1[0]=14.2;
    v1[1]=3.2;
    v1[2]=9.3;
    //v1[3]=0.0;
    
 
    v2[0]=1.2;
    v2[1]=32.2;
    v2[2]=12.0;
    //v2[3]=0.0;
    
    
    clock_t t;
    float time;
    type* res;
 
    t = clock();
    //sum(v1,v2,v3);
    euclidean_dist_sse(v1,v2, res);
    //type dist = euclidean_distance(v1,v2);
    
    t = clock() - t;
    time = ((float)t)/CLOCKS_PER_SEC;
 
       type* test = (type*) _mm_malloc(3 * sizeof(type), 16);
 
if (test == NULL) {
    printf("Memory allocation failed!\n");
} else {
    _mm_free(test);
}
 
    printf("time con assembly= %.3f secs\n", time);
    printf("DISTANCE: %f\n", *res);
}