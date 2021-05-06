#include "driver.h"	    // Include system files and define variables
#include "cFiles.h"
#include "limits.h"

struct indeces {
    int x;
    int y;
};

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

// Moser–De Bruijn sequence lookup and inverse lookup taken from:
// https://gist.github.com/JLChnToZ/ec41b1b45987d0e1b40ceabc13920559

int z_order_lookup(int x, int y) {
    if(x < 0 || y < 0) return 0;
    int result = 0;
    int mask;
    int offset = 0;
    
    for(mask = 1; x >= mask || y >= mask; mask <<= 1) {
        result |= (x & mask) << offset++ | (y & mask) << offset;
    }
    
    return result;
}

struct indeces z_order_inverse_lookup(int value) {
    struct indeces idx;
    idx.x = 0;
    idx.y = 0;
    int offset = 0;
    int mask;
    
    for(mask = 1; value == offset + 1 || value >= 1 << offset + 1; mask <<= 1) {
        idx.x |= (value >> offset++) & mask;
        idx.y |= (value >> offset) & mask;
    }
    return idx;
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
    
    int matmul = true;              /* Sequential Matrix Multiplication */
    int strassen = true;           /* Sequential Strassen Algorithm */
    int matmul_parallel = true;    /* Parallel Matrix Multiplication */
    int strassen_parallel = true;  /* Parallel Strassen Algorithm */
    
    double **A, **B, **C_seq, **C, **A_new, **B_new;
    double *ptrA, *ptrB, *ptrA_new, *ptrB_new;
    int i, j, run;
    
    int dim;  // dimension of matrix
    
    double mt1, mt2; // Timing variables
    float t_bs, t;
    
    sscanf(argv[1], "%d", &dim);
    
    t_bs = -1;
    
    printf("Allocating memory \n");
    
    A = allocate_matrix(dim);
    B = allocate_matrix(dim);
    
    srand(time(NULL));
    
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
            
            A_new = allocate_matrix(new_dim);
            B_new = allocate_matrix(new_dim);
            
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
    
    if (strassen) {
        for(run = 0; run < n_runs; run++) {
            //mt1 = omp_get_wtime();
            
            if (new_dim != dim) {
                C = sequential_strassen(A_new, B_new, new_dim, &t);
            } else {
                C = sequential_strassen(A, B, new_dim, &t);
            }
            
            //mt2 = omp_get_wtime();

            verify_matmul(C, C_seq, dim);
            
            //*** Capture best run
            
            if ((t_bs < 0) || (t < t_bs))
                t_bs = t;
        }
        
        if (dim <= 10) {
            printf("Strassen: \n");
            printf("C = ");
            printmatrix(C, dim);
        }
        
        printf("Done computing \n");
        printf("Strassen matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }
    
    t_bs = -1;

    
    if (matmul_parallel) {
        for(run = 0; run < n_runs; run++) {
            mt1 = omp_get_wtime();
            
            C = parallel_matmul(A, B, dim);
            
            mt2 = omp_get_wtime();
            
            //*** Capture best run
            
            if ((t_bs < 0) || (mt2 - mt1 < t_bs))
                t_bs = mt2 - mt1;
        }
        
        if (dim <= 10) {
            printf("Parallel: \n");
            printf("C = ");
            printmatrix(C_seq, dim);
        }
        
        printf("Done computing \n");
        printf("Parallel matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }

    t_bs = -1;
    
    if (strassen_parallel) {
        for(run = 0; run < n_runs; run++) {
            //mt1 = omp_get_wtime();
            
            if (new_dim != dim) {
                C = parallel_strassen(A_new, B_new, new_dim, &t);
            } else {
                C = parallel_strassen(A, B, new_dim, &t);
            }
            
            //mt2 = omp_get_wtime();

            verify_matmul(C, C_seq, dim);
            
            //*** Capture best run
            
            if ((t_bs < 0) || (t < t_bs))
                t_bs = t;
        }
        
        if (dim <= 10) {
            printf("Parallel Strassen: \n");
            printf("C = ");
            printmatrix(C, dim);
        }
        
        printf("Done computing \n");
        printf("Parallel Strassen matrix multiplication with %d x %d matrices took %f seconds\n", dim, dim, t_bs);
    }


    // TODO: fix memory deallocation
    
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
