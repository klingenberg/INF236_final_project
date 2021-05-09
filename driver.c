#include "driver.h"	    // Include system files and define variables
#include "cFiles.h"
#include "limits.h"

void printmatrix(double ** A, int n) {
    int i, j;
    for(i = 0; i < n; i++) {
        if (i == 0) printf("\t[");
        else printf("\t ");
        printf("[");
        for(j = 0; j < n; j++) {
            printf("%.4f", A[i][j]);
            if (j < n - 1) {
                printf(", ");
            }
        }
        printf("]");
        if (i == n - 1) printf("]");
        printf("\n");
    }
}

double ** allocate_matrix(int dim, double *ptr) {
    double **C;

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

int verify_matmul(double ** X, double **T, int dim) {

    int i, j;
    double eps = 0.000001;
    
    for(i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (fabs(X[i][j] - T[i][j]) > eps) {
                printf("Sequential: C[%d][%d] = %f, Parallel: C[%d][%d] = %f \n", i, j, T[i][j], i, j, X[i][j]);
                return false;
            }
        }
    }
          
    return true;
}

int main(int argc, char *argv[]) {
    
    int matmul = false;              /* Sequential Matrix Multiplication */
    int strassen = false;           /* Sequential Strassen Algorithm */
    int matmul_parallel = true;    /* Parallel Matrix Multiplication */
    int strassen_parallel_2_layers = true;  /* Parallel Strassen Algorithm */
    int strassen_parallel = true;  /* Parallel Strassen Algorithm */
    
    double **A, **B, **C_seq, **C, **A_new, **B_new, **C_new;
    double *ptrA, *ptrB, *ptrA_new, *ptrB_new, *ptrC_seq, *ptrC, *ptrC_new;
    int i, j, run;
    
    int dim;  // dimension of matrix
    
    double mt1, mt2; // Timing variables
    float t_bs, t;
    
    sscanf(argv[1], "%d", &dim);
    
    t_bs = -1;
    
    printf("Allocating memory \n");
    
    A = allocate_matrix(dim,ptrA);
    B = allocate_matrix(dim,ptrB);
    C = allocate_matrix(dim,ptrC);


    //srand(time(NULL));
    
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            A[i][j] = (double) rand() / (double) RAND_MAX;
            B[i][j] = (double) rand() / (double) RAND_MAX;
        }
    }
    
    if (dim <= 10) {
        printf("A = ");
        printmatrix(A, dim);
        printf("B = ");
        printmatrix(B, dim);
    }
    

    int new_dim = dim;
    if (strassen || strassen_parallel) {
        // check if dimension is a power of 2
        if ((dim & (dim - 1)) != 0) {
            new_dim = 1;
            while(new_dim < dim) new_dim *= 2;
            
            A_new = allocate_matrix(new_dim,ptrA_new);
            B_new = allocate_matrix(new_dim,ptrB_new);
            C_new = allocate_matrix(new_dim,ptrC_new);
            
            for (i = 0; i < dim; i++) {
                for (j = 0; j < dim; j++) {
                    A_new[i][j] = A[i][j];
                    B_new[i][j] = B[i][j];
                }
            }
        }
    }
    
    if (matmul) {
        for(run = 0; run < n_runs; run++) {
            mt1 = omp_get_wtime();
            
            C_seq = allocate_matrix(dim,ptrC_seq);
            sequential_matmul(C_seq, A, B, dim);
            
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
    
    if (strassen) {

        for(run = 0; run < n_runs; run++) {

            if (new_dim != dim) {
                sequential_strassen(C_new, A_new, B_new, new_dim, &t);
                //verify_matmul(C_new, C_seq, dim);
            } else {
                sequential_strassen(C, A, B, new_dim, &t);
                //verify_matmul(C, C_seq, dim);
            }

            
            
            //*** Capture best run
            
            if ((t_bs < 0) || (t < t_bs))
                t_bs = t;
        }
        
        /*
        if (dim <= 10) {
            printf("Strassen: \n");
            printf("C = ");
            printmatrix(C, dim);
        }*/
        
        printf("Done computing \n");
        printf("Strassen matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }
    
    t_bs = -1;

    
    if (matmul_parallel) {
        for(run = 0; run < n_runs; run++) {
            mt1 = omp_get_wtime();
            
            parallel_matmul(C, A, B, dim);
            
            mt2 = omp_get_wtime();
            
            //verify_matmul(C, C_seq, dim);
            
            //*** Capture best run
            
            if ((t_bs < 0) || (mt2 - mt1 < t_bs))
                t_bs = mt2 - mt1;
        }
        /*
        if (dim <= 10) {
            printf("Parallel: \n");
            printf("C = ");
            printmatrix(C_seq, dim);
        }
        */
        
        printf("Done computing \n");
        printf("Parallel matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }

    t_bs = -1;
    
    if (strassen_parallel_2_layers) {
        for(run = 0; run < n_runs; run++) {
            
            if (new_dim != dim) {
                parallel_strassen_2_layers(C_new, A_new, B_new, new_dim, &t);
                //verify_matmul(C_new, C_seq, dim);
            } else {
                parallel_strassen_2_layers(C, A, B, new_dim, &t);
                //verify_matmul(C, C_seq, dim);
            }
            
            //*** Capture best run
            
            if ((t_bs < 0) || (t < t_bs))
                t_bs = t;
        }
        /*
        if (dim <= 10) {
            printf("Parallel Strassen: \n");
            printf("C = ");
            printmatrix(C, dim);
        }
        */
        
        printf("Done computing \n");
        printf("Parallel 2-layers Strassen matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }

    t_bs = -1;
    
    if (strassen_parallel) {
        for(run = 0; run < n_runs; run++) {
            
            if (new_dim != dim) {
                parallel_strassen(C_new, A_new, B_new, new_dim, &t);
                //verify_matmul(C_new, C_seq, dim);
            } else {
                parallel_strassen(C, A, B, new_dim, &t);
                //verify_matmul(C, C_seq, dim);
            }
            
            //*** Capture best run
            
            if ((t_bs < 0) || (t < t_bs))
                t_bs = t;
        }
        /*
        if (dim <= 10) {
            printf("Parallel Strassen: \n");
            printf("C = ");
            printmatrix(C, dim);
        }
        */
        
        printf("Done computing \n");
        printf("Parallel Strassen matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }

    
    free(A);
    free(ptrA);
    free(B);
    free(ptrB);

    if(matmul) {
        free(C_seq);
        free(ptrC_seq);
    }

    free(C);
    free(ptrC);

    if (new_dim != dim) {
        free(A_new);
        free(B_new);
        free(C_new);
        free(ptrC_new);
    }

}
