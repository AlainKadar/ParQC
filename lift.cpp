/* Variable dictionary
 *
 * N - Number of points
 * R - Number of real dimensions
 * M - Number of higher dimensions
 * K - Number of neighbours * 2
 * J - Number of neighbours of single particle
 *
 * P.txt - Positions input file (flattened)
 * P - List of positions (flattened)
 * NL.txt - Neighbour list input file (flattened)
 * NL - Neighbour list (flattened)
 *
 * H - Hyperlattice matrix
 *
 * i - Only ever a for loop counter
 */

#include <omp.h>
#include <math.h>
#include <iostream>
#include <queue>
#include <stdio.h>

#include "printers.h"

//#define N 102891
#define R 3
#define M 5
//#define K 566188

#define MAX_NEIGHBOURS 5
#define MAX_ITER 500000

using namespace std;

//! Finds the location of a given argument in an array
/*! \param arr The array to search
    \param arg The argument to find in arr
    \param size The size of the array

    \returns The location of arg in arr.
*/
int ind_finder(int arr[], int arg, int size) {
    for (int i=0; i<size; i++){
        if (arr[i]==arg) { return i; }
    }
    throw; 
}

//! Finds maximum values of dot products in the dot product matrix
/*! \param dot_prod The matrix of dot products
    \param dot_max The array of maximum dot product values
    \param dot_max_ind The index locations of each of the maximum values
    \param J The number of neighbours (= number of rows in matrix)

    \returns void.
*/
void amax(float dot_prod[], float dot_max[], int dot_max_ind[], int J) {
    float _max;
    int _max_ind;
    int col = 0;
    int row = 0;
    for (int i=0; i<J*M; i++){
        if (col==0 || fabs(dot_prod[i]) > fabs(_max)) {
            _max=dot_prod[i];
            _max_ind=col;
        }

        if (col==(M-1)) {
            dot_max[row] = _max;
            dot_max_ind[row] = _max_ind;
            col=0;
            row++;
            continue;
        }
        col++;
    }
}

//! Matrix operation equivalent to A x B_transpose
/*! \param vij The matrix of positional difference vectors
    \param Q The basis set
    \param dot_prod The matrix of dot products, resulting from vij x Q_T
    \param J The number of neighbours (= number of rows in matrix)

    \returns void.
*/
void mat_mul(float vij[], float Q[], float dot_prod[], int J) {
    int n,m;
    for (int i=0; i<J*M; i++){
        dot_prod[i] = 0;
        n = floor(double(i)/double(M));
        m = i%M;
        for (int j=0; j<R; j++){
            dot_prod[i] = dot_prod[i] + Q[m*R+j]*vij[n*R+j];
        }
    }
}

//! Tests if element is in given array
/*! \param q The element to search for
    \param arr The array to search through
    \param arr_size The size of the array

    \returns true if q in arr. Otherwise returns false.
*/
bool in(int q, int arr[], int arr_size) {
    for (int i=0; i<arr_size; i++) {
        if (arr[i]==q) { return true; }
    }
    return false;
}

int main (int argc, char** argv) {
    

    char *_name = argv[1];
    string name = _name;
    
    char *A_string = argv[2];
    char *B_string = argv[3];
    char *C_string = argv[4];
    char *O_string = argv[5];
    char *S_string = argv[6];

    int A = atoi(A_string); //Some linters will label these as unused
    int B = atoi(B_string); //variables. However they are used by the omp
    int C = atoi(C_string); //directives.
    int O = atoi(O_string);
    float scale = atof(S_string);

    omp_set_nested(1);
    omp_set_max_active_levels(3);
    omp_lock_t Ref_lock, H_lock;
    omp_init_lock(&Ref_lock);
    omp_init_lock(&H_lock);

    float _temp_f;
    int _temp_i;

    string fnameP = "P_" + name + ".txt";
    string fnameNL = "NL_" + name + ".txt";
    string fnameH = "H_" + name + ".txt";
    int c;
    FILE *buffP = fopen(fnameP.c_str(), "r");
    FILE *buffNL = fopen(fnameNL.c_str(), "r");
    
    int i = 0;
    while (c != EOF) {
        fscanf(buffP, "%f", &_temp_f);
        c = getc(buffP);
        i++;
    }
    int N = int(i/R);

    
    int d;
    i = 0;
    while (d != EOF) {
        fscanf(buffNL, "%i", &_temp_i);
        d = getc(buffNL);
        i++;
    }
    int K = i;
    
    float P[(N+10)*R];
    int NL[K+10];
    int H[(N+10)*M] = { 0 };
    rewind(buffP);
    rewind(buffNL);
     
    for (int i=0; i<N*R; i++) {
        fscanf(buffP, "%f", &P[i]);
    }
    
    for (int i=0; i<K; i++) {
        fscanf(buffNL, "%i", &NL[i]);
    }

    fclose(buffP);
    fclose(buffNL);
    float PI=3.1415;
    //float scale=0.917; //For ideal_tiling.gsd
    //scale=0.964; //For flattened_simulation.gsd
   
    float Q[M*R] = {0, 1*scale, 0,
               sin(2*PI/5)*scale, cos(2*PI/5)*scale, 0,
               sin(4*PI/5)*scale, cos(4*PI/5)*scale, 0,
               sin(6*PI/5)*scale, cos(6*PI/5)*scale, 0,
               sin(8*PI/5)*scale, cos(8*PI/5)*scale, 0};

    int origin_ind = O;
    int n = 0;
    int k = 1;
    int h = M;
    int visited[N] = { N + 1 }; //Initialise such that the array contains a value that refers to no particle in P
    int H_ind[N];
    H_ind[0] = origin_ind;
    visited[0] = origin_ind;
    int visited_size = 1;
    queue<int> Refs;
    Refs.push(origin_ind);
   
    bool finished=false;
    int ref_H[M] = { 0, 0, 0, 0, 0 };
    int ref_H_ind;
    int J;
    int ind_ref;
    int nl_mod_ind [MAX_NEIGHBOURS*2]; //Initialise such that it could contain  maximum number of neighbours.

    double start, end;
    start = omp_get_wtime();
    while (!finished)
    {
        #pragma omp parallel for private(ind_ref, ref_H, ref_H_ind, J, nl_mod_ind) num_threads(A)
        for (int m=0; m<MAX_ITER; m++) 
        {
            if (visited_size>=N) {
                m=MAX_ITER;
                finished=true;}
            omp_set_lock(&Ref_lock);
            if (Refs.empty()) {
                omp_unset_lock(&Ref_lock);
                continue;
            } else

            {
                ind_ref = Refs.front();
                Refs.pop();

                ref_H_ind = ind_finder(H_ind, ind_ref, k);
                ref_H[0] = H[M*ref_H_ind];
                ref_H[1] = H[M*ref_H_ind+1];
                ref_H[2] = H[M*ref_H_ind+2];
                ref_H[3] = H[M*ref_H_ind+3];
                ref_H[4] = H[M*ref_H_ind+4];
                
                int ind_size; // = len(ind)/R
               
                int temp=0;
                #pragma omp parallel num_threads(B)
                {
                    int i, id, nthrds;
                    id = omp_get_thread_num();
                    nthrds = omp_get_num_threads();
                    for (i=2*id; i<K; i=i+2*nthrds) {
                        if (NL[i] == ind_ref) {
                            //if (!in(NL[i+1], visited, visited_size)){
                            //bool in(int q, int arr[], int arr_size) {
                            if (visited_size<=2000) {
                                if (!in(NL[i+1], visited, visited_size)){
                                    #pragma omp critical
                                    {
                                        nl_mod_ind[temp] = i;
                                        nl_mod_ind[temp+1]=i+1;
                                        temp = temp+2;
                                    }
                                }
                            } else {
                                bool found=false;
                                #pragma omp parallel num_threads(C)
                                {
                                    int id_int, nthrds_int;
                                    id_int = omp_get_thread_num();
                                    nthrds_int = omp_get_num_threads();
                                    for (int i_int=id_int; i_int<visited_size; i_int=i_int+nthrds_int) {
                                        if (NL[i+1]==visited[i_int]) {
                                            found=true;
                                        }
                                    }
                                }
                                if (!found) {
                                    #pragma omp critical
                                    {
                                        nl_mod_ind[temp] = i;
                                        nl_mod_ind[temp+1]=i+1;
                                        temp = temp+2;
                                    }
                                }
                            }
                        }
                    }
                    ind_size = temp;
                }
                if (ind_size==0) {     
                    omp_unset_lock(&Ref_lock);
                    continue;
                } else 
                {
                    J = int(ind_size/2);
                    for (int neighbour=0; neighbour<J; neighbour++) {
                        visited[visited_size] = NL[nl_mod_ind[2*(neighbour)+1]];
                        visited_size++;
                                        }
                    omp_unset_lock(&Ref_lock);
                       
                    float vij[J*R];
                    for (int i=0; i<J; i++){
                        vij[i*R]   = P[R*ind_ref]   - P[R*NL[nl_mod_ind[2*i+1]]];
                        vij[i*R+1] = P[R*ind_ref+1] - P[R*NL[nl_mod_ind[2*i+1]]+1];
                        vij[i*R+2] = P[R*ind_ref+2] - P[R*NL[nl_mod_ind[2*i+1]]+2];
                    }
                    float dot_prod[J*M];
                    mat_mul(vij, Q, dot_prod, J);

                    float dot_max[J];
                    int dot_max_ind[J];
                    amax(dot_prod, dot_max, dot_max_ind, J);
                    

                    omp_set_lock(&H_lock);
                    for (int neighbour=1; neighbour<J+1; neighbour++) {
                        for (int column=0; column<M; column++) {
                            if (column==dot_max_ind[neighbour-1]) {
                                if (dot_max[neighbour-1] > 0) { H[h] = ref_H[column] + 1; }
                                else { H[h] = ref_H[column] -1; }
                            } else {
                                H[h] = ref_H[column];
                            }
                            h++;
                        }

                        
                        H_ind[k] = NL[nl_mod_ind[2*(neighbour-1)+1]];
                        k++;
                        

                        Refs.push(NL[nl_mod_ind[2*(neighbour-1)+1]]);
                    }
                    omp_unset_lock(&H_lock);
                    n++;
                }
            }
        }
        finished=true;
    }

    end = omp_get_wtime();
    FILE *buffH = fopen(fnameH.c_str(), "w");
    for (int i=0; i<M*k; i++) {
        fprintf(buffH, "%i \n", H[i]);
    }

    fclose(buffH);
    omp_destroy_lock(&Ref_lock);
    omp_destroy_lock(&H_lock);
    printf("Time taken was %f\n", end-start);
}



