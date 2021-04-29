#include "driver.h"	    // Include system files and define variables
#include "cFiles.h"

void printmatrix(int ** A, int n) {
    int i, j;
    for(i = 0; i < n; i++) {
        if (i == 0) printf("\t[");
        else printf("\t ");
        printf("[");
        for(j = 0; j < n; j++) {
            printf("%d", A[i][j]);
            if (j < n - 1) {
                printf(", ");
            }
        }
        printf("]");
        if (i == n - 1) printf("]");
        printf("\n");
    }
}

int ** allocate_matrix(int dim) {
    int **C;
    int * ptr;
    ptr = (int*) malloc(sizeof(int) * dim * dim);
    C = (int**) malloc(sizeof(int *) * dim);
    
    int i;
    for(i = 0; i < dim; i++) C[i] = ptr + i * dim;
    
    if (C == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }
    return C;
}

int verify_matmul(int ** X, int **T, int dim) {

    int i, j;
    
    for(i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (X[i][j] != T[i][j]) {
                printf("Sequential: C[%d][%d] = %d, Parallel: C[%d][%d] = %d \n", i, j, T[i][j], i, j, X[i][j]);
                return false;
            }
        }
    }
          
    return true;
}

int main(int argc, char *argv[]) {
    
    int matmul = true;              /* Sequential Matrix Multiplication */
    int strassen = true;           /* Sequential Strassen Algorithm */
    int matmul_parallel = false;    /* Parallel Matrix Multiplication */
    int strassen_parallel = false;  /* Parallel Strassen Algorithm */
    
    int **A, **B, **C_seq, **C, **A_new, **B_new;
    int i, j, run;
    
    int dim;  // dimension of matrix
    
    double mt1, mt2; // Timing variables
    float t_bs;
    
    printf("Give matrix dimension \n");
    scanf("%d", &dim);
    
    t_bs = -1;
    
    printf("Allocating memory \n");
    
    A = allocate_matrix(dim);
    B = allocate_matrix(dim);
    
    srand(time(NULL));
    
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            A[i][j] = rand() % 20;
            B[i][j] = rand() % 20;
        }
    }
    
    if (dim <= 10) {
        printf("A = ");
        printmatrix(A, dim);
        printf("B = ");
        printmatrix(B, dim);
    }
    
    if (matmul) {
        for(run = 0; run < n_runs; run++) {
            mt1 = omp_get_wtime();
            
            C_seq = sequential_matmul(A, B, dim);
            
            mt2 = omp_get_wtime();
            
            //*** Capture best run
            
            if ((t_bs < 0) || (mt2 - mt1 < t_bs))
                t_bs = mt2 - mt1;
        }
        
        if (dim <= 10) {
            printf("Sequential: \n");
            printf("C = ");
            printmatrix(C_seq, dim);
        }
        
        printf("Done computing \n");
        printf("Simple matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }
    
    t_bs = -1;
    
    int new_dim = dim;
    if (strassen) {
        // check if dimension is a power of 2
        if ((dim & (dim - 1)) != 0) {
            new_dim = 1;
            while(new_dim < dim) new_dim *= 2;
            
            A_new = allocate_matrix(new_dim);
            B_new = allocate_matrix(new_dim);
            
            for (i = 0; i < dim; i++) {
                for (j = 0; j < dim; j++) {
                    A_new[i][j] = A[i][j];
                    B_new[i][j] = B[i][j];
                }
            }
        }
        
        for(run = 0; run < n_runs; run++) {
            mt1 = omp_get_wtime();
            
            if (new_dim != dim) {
                C = allocate_matrix(new_dim);
                sequential_strassen(C, A_new, B_new, new_dim);
            } else {
                C = allocate_matrix(dim);
                sequential_strassen(C, A, B, new_dim);
            }
            
            mt2 = omp_get_wtime();

            verify_matmul(C, C_seq, dim);
            
            //*** Capture best run
            
            if ((t_bs < 0) || (mt2 - mt1 < t_bs))
                t_bs = mt2 - mt1;
        }
        
        if (dim <= 10) {
            printf("Strassen: \n");
            printf("C = ");
            printmatrix(C, dim);
        }
        
        printf("Done computing \n");
        printf("Strassen matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }
    
    free(A);
    free(B);
    free(C_seq);

    if (strassen) {
        free(C);
    }

    if (new_dim != dim) {
        free(A_new);
        free(B_new);
    }

}
