#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define n_runs 10

void printmatrix(double ** A, int m, int n) {
    
    int i, j;
    for(i = 0; i < m; i++) {
        if (i == 0) printf("\t[");
        else printf("\t ");
        printf("[");
        for(j = 0; j < n; j++) {
            printf("%.4f", A[i][j]);
            if (j < m - 1) {
                printf(", ");
            }
        }
        printf("]");
        if (i == m - 1) printf("]");
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    
    double **A;
    double **B;
    double **C;
    double *ptr1,*ptr2,*ptr3;
    int i,j,k,l;
    
    int dim;  // dimension of matrix
    
    double mt1,mt2; // Timing variables
    float t_bs;
    
    printf("Give matrix dimension \n");
    scanf("%d",&dim);
    
    printf("Allocating memory \n");
    
    A = malloc(sizeof(double *)*dim);
    ptr1 = malloc(sizeof(double)*dim*dim);
    B = malloc(sizeof(double *)*dim);
    ptr2 = malloc(sizeof(double)*dim*dim);
    C = malloc(sizeof(double *)*dim);
    ptr3 = malloc(sizeof(double)*dim*dim);
    
    if (C == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }
    
    for(i=0;i<dim;i++) {
        A[i] = ptr1 + (i*dim);
        B[i] = ptr2 + (i*dim);
        C[i] = ptr3 + (i*dim);
    }
    
    t_bs = -1.0;
    
    srand(time( NULL));
    
    for(i = 0; i < dim; i++) {
        for(j = 0;j < dim; j++) {
            A[i][j] = (double) rand() / (double) RAND_MAX;
            B[i][j] = (double) rand() / (double) RAND_MAX;;
        }
    }
    
    if (dim <= 10) {
        printf("A = ");
        printmatrix(A, dim, dim);
        printf("B = ");
        printmatrix(B, dim, dim);
    }
    
    printf("Starting the Matrix multiplication algorithm \n");
    
    for(l = 0; l < n_runs; l++) {
        mt1 = omp_get_wtime();
        
        // *********************************
        // Sequential matrix multiplication
        // *********************************
        
        for(i = 0;i < dim; i++) {
            for(j = 0;j < dim; j++) {
                C[i][j] = 0.0;
                for(k = 0;k < dim; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                } // k
            } // j
        } // i
        
        
        mt2 = omp_get_wtime();
        
        //*** Capture best run
        
        if ((t_bs < 0) || (mt2-mt1 < t_bs))
            t_bs = mt2-mt1;
    }
    
    if (dim <= 10) {
        printf("C = ");
        printmatrix(C, dim, dim);
    }
    
    printf("Done computing \n");
    printf("Matrix multiplication with %d x %d matrices took %f seconds\n",dim,dim,t_bs);
}

