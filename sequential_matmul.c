double ** sequential_matmul(double ** A, double ** B, int dim) {
    
    int i, j, k;
    double **C = allocate_matrix(dim);
    
    // *********************************
    // Sequential Matrix Multiplication
    // *********************************
    
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            C[i][j] = 0.0;
        } // j
        for(k = 0; k < dim; k++) {
            for (j = 0; j < dim; j++) {
                C[i][j] += A[i][k] * B[k][j];
            } // j
        } // k
    } // i
    
    return C;
}
