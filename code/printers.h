#include <iostream>


void print_float_matrix(float A[], int len, int _size){
    for (int i = 0; i < _size; i++){
        if (i%len==(len-1)){
            printf("%f, \n", A[i]);
        } else {
            printf("%f, ", A[i]);
        }
    }
}

void print_int_matrix(int A[], int len, int _size){
    for (int i = 0; i < _size; i++){
        if (i%len==(len-1)){
            printf("%i, \n", A[i]);
        } else {
            printf("%i, ", A[i]);
        }
    }
}

void print_float_array(float p[], int _size){
    for (int i=0; i<_size; i++){
        printf("%f, ", p[i]);
    }
    printf("\n");
}
void print_int_array(int p[], int _size){
    for (int i=0; i<_size; i++){
        printf("%i, ", p[i]);
    }
    printf("\n");
}

