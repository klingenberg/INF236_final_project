#include "driver.h"	    // Include system files and define variables
#include "cFiles.h"

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
    
    int matmul = true;              /* Sequential Matrix Multiplication */
    int strassen = true;           /* Sequential Strassen Algorithm */
    int matmul_parallel = false;    /* Parallel Matrix Multiplication */
    int strassen_parallel = false;  /* Parallel Strassen Algorithm */
    
    double **A, **B, **C;
    double *ptrA,*ptrB;
    int i, j, run;
    
    int dim;  // dimension of matrix
    
    double mt1, mt2; // Timing variables
    float t_bs;
    
    printf("Give matrix dimension \n");
    scanf("%d",&dim);
    
    printf("Allocating memory \n");
    
    A = malloc(sizeof(double *)*dim);
    ptrA = malloc(sizeof(double)*dim*dim);
    B = malloc(sizeof(double *)*dim);
    ptrB = malloc(sizeof(double)*dim*dim);
    
    for(i = 0; i < dim; i++) {
        A[i] = ptrA + (i*dim);
        B[i] = ptrB + (i*dim);
    }
    
    srand(time(NULL));
    
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
    
    if (matmul) {
        for(run = 0; run < n_runs; run++) {
            mt1 = omp_get_wtime();
            
            C = sequential_matmul(A, B, dim);
            
            mt2 = omp_get_wtime();
            
            //*** Capture best run
            
            if ((t_bs < 0) || (mt2-mt1 < t_bs))
                t_bs = mt2-mt1;
        }
    }
    
    if (dim <= 10) {
        printf("C = ");
        printmatrix(C, dim, dim);
    }
    
    if (strassen) {
        for(run = 0; run < n_runs; run++) {
            mt1 = omp_get_wtime();
            
            C = sequential_strassen(A, B, dim);
            
            mt2 = omp_get_wtime();
            
            //*** Capture best run
            
            if ((t_bs < 0) || (mt2-mt1 < t_bs))
                t_bs = mt2-mt1;
        }
    }
    
    printf("Done computing \n");
    printf("Matrix multiplication with %d x %d matrices took %f seconds\n",dim,dim,t_bs);
}
