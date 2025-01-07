#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <malloc.h>
#include <stdint.h>
 
#define type        double
#define MATRIX      type*
#define REAL_MATRIX type**
#define VECTOR      type*
 
#define random() (((type) rand())/RAND_MAX)
 
type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};      // hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};      // volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};       // charge
 
typedef struct {
    char* seq;      // sequenza di amminoacidi
    int N;          // lunghezza sequenza
    unsigned int sd;    // seed per la generazione casuale
    type to;        // temperatura INIZIALE
    type alpha;     // tasso di raffredamento
    type k;     // costante
    VECTOR hydrophobicity; // hydrophobicity
    VECTOR volume;      // volume
    VECTOR charge;      // charge
    VECTOR phi;     // vettore angoli phi
    VECTOR psi;     // vettore angoli psi
    type e;     // energy
    int display;
    int silent;
 
} params;
 
 
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
 

void gen_rnd_mat(VECTOR v, int N){
    int i;
 
    for(i=0; i<N; i++){
        // Campionamento del valore + scalatura
        v[i] = (random()*2 * M_PI) - M_PI;
    }
}
 
void sub(VECTOR v1, VECTOR v2, VECTOR res){
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
 
void sum(VECTOR v1, VECTOR v2, VECTOR res){
    for(int i=0;i< 4; i++){
        res[i]=v1[i]+v2[i];
    }
}
 

 
void prodotto_vettore_scalare(VECTOR v, type s){
    for(int i=0;i<4;i++){
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
 
 
void prodotto_vettore_matrice(VECTOR v, REAL_MATRIX m, VECTOR res){
    // Esegui il prodotto tra il vettore e la matrice
    for (int i = 0; i < 3; i++) { // Itera sulle colonne della matrice
        for (int j = 0; j < 3; j++) { // Itera sulle righe della matrice
            res[i] += v[j] * m[j][i];
        }
    }
}
 
void mul(VECTOR v1, VECTOR v2){
    for(int i=0;i<4;i++){
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
    type norm = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2] + axis[3]*axis[3];
    norm = sqrtf(norm); // Calcolo della radice quadrata
    // Normalizza il vettore axis
    axis[0] /= norm;
    axis[1] /= norm;
    axis[2] /= norm;
    axis[3] /= norm;
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
 
 
void rotation(VECTOR axis, type theta, REAL_MATRIX rotation_matrix){
    type res;
    type a, b, c, d;
    prodotto_scalare(axis, axis, &res);
 
    div_scalareC(axis, res);
    prodotto_vettore_scalare(axis, approx_sin(theta / 2.0f));
    prodotto_vettore_scalare(axis, -1);
 
    a=approx_cos(theta / 2.0f);
    b=axis[0];
    c=axis[1];
    d=axis[2];
 
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
}
 
 
MATRIX backbone(char* s, int n, VECTOR phi, VECTOR psi){
 
    //distanze atomi in Angstrom
    type r_ca_n = 1.46;
    type r_ca_c = 1.52;
    type r_c_n = 1.33;
    
    //angoli atomi in radianti
    type theta_ca_c_n = 2.028;
    type theta_c_n_ca = 2.124;
    type theta_n_ca_c = 1.940;
    
    //allocazione matrice coords
    MATRIX coords = alloc_matrix( n*3, 3);
    
    int i, j;
    
    //inizializzazione coordinate N primo amminoacido
    coords[0]=0;
    coords[1]=0;
    coords[2]=0;
 
    //iniziaizzazione coordinate C_alpha primo amminoacido
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
        
        int idx=i*3*3; //calcolo indice base amminoacido
        
        if (i>0) {
 
            //posiziona N
        
            //popolo i vettori
            coords_c_alpha[0] = coords[(i-1)*3*3+3];
            coords_c_alpha[1] = coords[(i-1)*3*3+4];
            coords_c_alpha[2] = coords[(i-1)*3*3+5];
            coords_c_alpha[3] = 0;
 
            coords_c[0] = coords[(i-1)*3*3+6];
            coords_c[1] = coords[(i-1)*3*3+7];
            coords_c[2] = coords[(i-1)*3*3+8];
            coords_c[3] = 0;
            
            // C-C_alpha
            //sub(coords_c, coords_c_alpha, v1);
            v1[0] = coords_c[0]-coords_c_alpha[0];
            v1[1] = coords_c[1]-coords_c_alpha[1];
            v1[2] = coords_c[2]-coords_c_alpha[2];
            v1[3] = coords_c[3]-coords_c_alpha[3];
            
            //v1/||v1||
            //normalize_axis(v1);
            norm = sqrtf(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] + v1[3]*v1[3]);
            v1[0] = v1[0]/norm;
            v1[1] = v1[1]/norm;
            v1[2] = v1[2]/norm;
            v1[3] = v1[3]/norm;
 
            //costruzione matrice rotation
            
            rotation(v1, theta_c_n_ca, rotation_matrix);
            
            
            v[0]=0;
            v[1]=r_c_n;
            v[2]=0;
            
            newv[0]=0;
            newv[1]=0;
            newv[2]=0;
            newv[3]=0;
 
            prodotto_vettore_matrice(v, rotation_matrix, newv);
            
            //N = C + newv
            //sum(coords_c, newv, coords_n);
            coords_n[0] = coords_c[0]+newv[0];
            coords_n[1] = coords_c[1]+newv[1];
            coords_n[2] = coords_c[2]+newv[2];
            coords_n[3] = coords_c[3]+newv[3];
            
            //inserisco in coords le coordinate di N
            coords[idx] = coords_n[0];
            coords[idx+1] = coords_n[1];
            coords[idx+2] = coords_n[2];
 
            //posiziona C_alpha
            //VECTOR v2 = alloc_matrix(1,4);
 
            //N - C
            //sub(coords_n, coords_c, v2);
            v2[0] = coords_n[0]-coords_c[0];
            v2[1] = coords_n[1]-coords_c[1];
            v2[2] = coords_n[2]-coords_c[2];
            v2[3] = coords_n[3]-coords_c[3];
            
            //v2/||v2||
            //normalize_axis(v2);
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
 
            prodotto_vettore_matrice(v_, rotation_matrix, newv);
            
            //C_alpha = N + newv
            //sum(coords_n, newv, coords_c_alpha);
            coords_c_alpha[0] = coords_n[0]+newv[0];
            coords_c_alpha[1] = coords_n[1]+newv[1];
            coords_c_alpha[2] = coords_n[2]+newv[2];
            coords_c_alpha[3] = coords_n[3]+newv[3];
 
            //inserisco in coords le coordinate di C_alpha
            coords[idx+3] = coords_c_alpha[0];
            coords[idx+4] = coords_c_alpha[1];
            coords[idx+5] = coords_c_alpha[2];
        }
        //Posiziona C
    
        coords_n[0] = coords[idx];
        coords_n[1] = coords[idx+1];
        coords_n[2] = coords[idx+2];
        coords_n[3] = 0;
 
        coords_c_alpha[0] = coords[idx+3];
        coords_c_alpha[1] = coords[idx+4];
        coords_c_alpha[2] = coords[idx+5];
        coords_c_alpha[3] = 0;
    
        //C_alpha - N
        //sub(coords_c_alpha, coords_n, v3);
        v3[0] = coords_c_alpha[0] - coords_n[0];
        v3[1] = coords_c_alpha[1] - coords_n[1];
        v3[2] = coords_c_alpha[2] - coords_n[2];
        v3[3] = coords_c_alpha[3] - coords_n[3];
    
        //v3/||v3||
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
        
        prodotto_vettore_matrice(v, rotation_matrix, newv);
    
        //C = C_alpha + newv
        coords_c[0] = coords_c_alpha[0] +newv[0];
        coords_c[1] = coords_c_alpha[1] +newv[1];
        coords_c[2] = coords_c_alpha[2] +newv[2];
        coords_c[3] = coords_c_alpha[3] +newv[3];
        //inserisco in coords le coordinate di C
    
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
    for(i=0; i<n; i++){
        type alpha_dist = sqrtf((phi[i]-alpha_phi)*(phi[i]-alpha_phi)+((psi[i]-alpha_psi)*(psi[i]-alpha_psi)));
        type beta_dist = sqrtf((phi[i]-beta_phi)*(phi[i]-beta_phi)+(psi[i]-beta_psi)*(psi[i]-beta_psi));
        type min = fmin(alpha_dist, beta_dist);
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
 
    //estrapolo coordinate C_alpha
    for(i=0; i<n ; i++){
        
        int idx_i = i*3*3+3; //indirizzo base generico C_alpha_i
    
        coords_c_alpha_i[0] = coords[idx_i];
        coords_c_alpha_i[1] = coords[idx_i+1];
        coords_c_alpha_i[2] = coords[idx_i+2];
        coords_c_alpha_i[3] = 0;
 
        
        
        for(j=i+1; j<n ; j++){
 
            int idx_j=j*3*3+3; //indirizzo base generico C_alpha_j
 
 
            coords_c_alpha_j[0] = coords[idx_j];
            coords_c_alpha_j[1] = coords[idx_j+1];
            coords_c_alpha_j[2] = coords[idx_j+2];
            coords_c_alpha_j[3] = 0;
 
            type dist = 0;
            euclidean_dist(coords_c_alpha_i, coords_c_alpha_j, &dist);
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
 
        //estrapolo coordinate C_alpha
        for(i=0; i<n ; i++){
 
            int idx_i = i*3*3+3; //indirizzo base generico C_alpha_i
    
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
 
                euclidean_dist(coords_c_alpha_i, coords_c_alpha_j, &dist);
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
 
 
        //estrapolo coordinate C_alpha_i e C_alpha_j
        for(i=0; i<n ; i++){
            type density = 0;
            int pos_i = s[i]-65;
 
            int idx_i = i*3*3+3; //indirizzo base generico C_alpha_i
                    
            coords_c_alpha_i[0] = coords[idx_i];
            coords_c_alpha_i[1] = coords[idx_i+1];
            coords_c_alpha_i[2] = coords[idx_i+2];
            coords_c_alpha_i[3] = 0;
                
            for(j=0; j<n ; j++){
 
                int idx_j=j*3*3+3; //indirizzo base generico C_alpha_j
        
                
                coords_c_alpha_j[0] = coords[idx_j];
                coords_c_alpha_j[1] = coords[idx_j+1];
                coords_c_alpha_j[2] = coords[idx_j+2];
                coords_c_alpha_j[3] = 0;
                       
                type dist = 0;
                euclidean_dist(coords_c_alpha_i, coords_c_alpha_j, &dist);
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
 
        //calcolo delle componenti energetiche
        type rama = rama_energy(phi, psi);
        type hydro = hydrophobic_energy(s,n, coords);
        type elec = electrostatic_energy(s,n, coords);
        type pack = packing_energy(s,n, coords);
 
    VECTOR v1 = alloc_matrix(1,4);
    //popolo vettore componenti energetiche
    v1[0] = rama; v1[1] = hydro; v1[2] = elec; v1[3] = pack;
 
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
 
    //printf("rama:%f, hydro:%f, elec:%f, pack:%f\n totale:%f\n", rama, hydro, elec, pack, total);
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
    //exit(0);
 
    while(T > 0){
 
        //genera un vicino della soluzione corrente
        unsigned int seed = input->sd;
        int i = (int) (random()*input->N); //prova anche altra formula col seed!
        //calcola variazioni casuali
        type theta_phi = (random()*2*M_PI)-M_PI;
        type theta_psi = (random()*2*M_PI)-M_PI;
 
        //aggiorna i valori degli angoli
        input->phi[i] = input->phi[i] + theta_phi;
        input->psi[i] = input->psi[i] + theta_psi;
 
        //calcola la variazione dell'energia
        type new_E = energy(input->seq,input->N, input->phi, input->psi);
        type delta_E = new_E - E;
        //printf("energy = %f\n", new_E);
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
        T = input->to - sqrt((input->alpha)*t);
        /*printf("energy <%d>: %f\n ",t, E );
        if(t==10) exit(0);*/
        
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
    sprintf(fname_phi, "out64_%d_%d_phi.ds2", input->N, input->sd);
    save_out(fname_phi, input->phi, input->N);
    sprintf(fname_psi, "out64_%d_%d_psi.ds2", input->N, input->sd);
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
 
 