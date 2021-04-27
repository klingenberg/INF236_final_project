double ** sequential_matmul(double ** A, double ** B, int dim) {
    
    int i, j, k;
    double **C;
    double * ptr;
    C = malloc(sizeof(double *) * dim);
    ptr = malloc(sizeof(double) * dim * dim);
    
    for(i = 0; i < dim; i++) {
        C[i] = ptr + (i*dim);
    }
    
    if (C == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }
    
    // *********************************
    // Sequential matrix multiplication
    // *********************************
    
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            C[i][j] = 0.0;
            for(k = 0; k < dim; k++) {
                C[i][j] += A[i][k] * B[k][j];
            } // k
        } // j
    } // i
    
    return C;
}
