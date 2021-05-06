#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define true  1
#define false 0
#define n_runs 2

double * allocate_array(int n) {
    double *T;
    T = (double*) malloc(sizeof(double) * n);

    if (T == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }

    return T;
}

void add(double *C, double *A, double *B, int n) {
    int i;
    for(i = 0; i < n; i++) C[i] = A[i] + B[i];
}

void parallel_add(double *C, double *A, double *B, int n) {
    int i;
    #pragma omp parallel for
    for(i = 0; i < n; i++) C[i] = A[i] + B[i];
}

int main(int argc, char *argv[]) {
    
    double *A, *B, *C;
    int i, j, run;
    
    int n;  // dimension of matrix
    
    double mt1, mt2; // Timing variables
    float t_bs, t;
    
    sscanf(argv[1], "%d", &n);

    int dim = n*n;
    
    t_bs = -1;
    
    A = allocate_array(dim);
    B = allocate_array(dim);
    
    srand(time(NULL));
    
    for(i = 0; i < dim; i++) {
            A[i] = (double) rand() / (double) RAND_MAX;
            B[i] = (double) rand() / (double) RAND_MAX;
    }
    
    C = allocate_array(dim);

    t_bs = -1;
    for(run = 0; run < n_runs; run++) {
        mt1 = omp_get_wtime();
        
        add(C, A, B, dim);
        
        mt2 = omp_get_wtime();
        
        //*** Capture best run
        
        if ((t_bs < 0) || (mt2 - mt1 < t_bs))
            t_bs = mt2 - mt1;
    }
    printf("Simple sum with %d elements took %f seconds\n", dim,  t_bs);

    t_bs = -1;
    for(run = 0; run < n_runs; run++) {
        mt1 = omp_get_wtime();
        
        parallel_add(C, A, B, dim);
        
        mt2 = omp_get_wtime();
        
        //*** Capture best run
        
        if ((t_bs < 0) || (mt2 - mt1 < t_bs))
            t_bs = mt2 - mt1;
    }
    printf("Parallel sum with %d elements took %f seconds\n", dim,  t_bs);
}