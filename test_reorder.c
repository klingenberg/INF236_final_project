#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define true  1
#define false 0

#define n_runs 1

double ** allocate_matrix(int dim) {
    double **C;
    double * ptr;
    ptr = (double*) malloc(sizeof(double) * dim * dim);
    C = (double**) malloc(sizeof(double *) * dim);
    
    int i;
    for(i = 0; i < dim; i++) C[i] = ptr + i * dim;
    
    if (C == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }
    return C;
}

double * allocate_array(int n) {
    double *T;
    T = (double*) malloc(sizeof(double) * n);

    if (T == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }

    return T;
}

void printmatrix(double ** A, int n) {
    int i, j;
    for(i = 0; i < n; i++) {
        if (i == 0) printf("\t[");
        else printf("\t ");
        printf("[");
        for(j = 0; j < n; j++) {
            printf("%.0f", A[i][j]);
            if (j < n - 1) {
                printf(", ");
            }
        }
        printf("]");
        if (i == n - 1) printf("]");
        printf("\n");
    }
}

void swap(double *xp,  double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void unshuffle(double *T, double *N, int n) {
    int i,j,k,l;

    for (i=0;i<n/4;i++) { 
        for (j=0;j<n;j++) {
            for (k=0;k<2;k++) {
                
                int row = 2*i;
                if (j>=n/2) row++;

                int col = j % (n/2) + k * (n/2);

                N[row*n+col] = T[i*n+k*(n*n/4)+j];
            }
        }
        
    }

    for (i = 0; i < n*n/2; i++) { 
        T[i] = N[i];
    }
}

void shuffle(double *T, double *N, int n) {
    int i,j,k,l;

    for (i=0;i<n/4;i++) { 
        for (j=0;j<n;j++) {
            for (k=0;k<2;k++) {
                
                int row = 2*i;
                if (j>=n/2) row++;

                int col = j % (n/2) + k * (n/2);

                N[i*n+k*(n*n/4)+j] = T[row*n+col];
            }
        }
        
    }

    for (i = 0; i < n*n/2; i++) { 
        T[i] = N[i];
    }
}

void mortonOrder(double *A, double *H, int n, int b) { 
    int size,submatrix;
    for(size = b; size <n;size = size *2) {
        int sub_size = (size*size);
        for(submatrix = 0; submatrix< n*n;submatrix = submatrix+sub_size*2) {
            unshuffle(A+submatrix,H,size*2);
        }
    }
}

void mortonOrderBack(double *A, double *H, int n, int b) {            
    int size,submatrix;
    for(size = n/2; size >=b;size = size /2) {
        int sub_size = (size*size);
        for(submatrix = 0; submatrix< n*n;submatrix = submatrix+sub_size*2) {
            shuffle(A+submatrix,H,size*2);
        }
    }
}

double * reorder_to_z_array(double ** X, int n, int min_square) {
    double *T;
    int dim = n*n;
    int i,j;
    T = allocate_array(dim);

    for(i=0;i<n;i++) for(j=0;j<n;j++) T[i*n+j] = X[i][j];

    double *help;
    help = allocate_array(dim);

    mortonOrder(T,help,n,min_square);

    return T;
}

void reorder_back_z_array(double ** X, double *T, int n, int min_square) {

    int dim = n*n;
    int i,j;

    double *help;
    help = allocate_array(dim);

    mortonOrderBack(T,help,n,min_square);

    for(i=0;i<n;i++) for(j=0;j<n;j++) X[i][j] = T[i*n+j];
}

int main(int argc, char *argv[]) {
    
    double **M;
    int i, j;
    
    int dim = 16; 
    int d = 4;

    M = allocate_matrix(dim);

    printf("Matrix: \n");
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            M[i][j] = i*dim+j;
            printf("%.0f ",M[i][j]);
        }printf("\n");
    }

    double *R;
    R = reorder_to_z_array(M, dim, d);


    printf("\nMorton:");
    for (j = 0; j < dim*dim; j++) {
        if (j % dim == 0) printf("\n");
        printf("%.0f ", R[j]);
    }printf("\n");

    reorder_back_z_array(M, R, dim, d);

    printf("\n");
    printf("Reordered back:\n");
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            M[i][j] = i*dim+j;
            printf("%.0f ",M[i][j]);
        }printf("\n");
    }

}