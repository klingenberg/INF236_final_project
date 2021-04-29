int ** sequential_matmul(int ** A, int ** B, int dim) {
    
    int i, j, k;
    int **C = allocate_matrix(dim);
    
    // *********************************
    // Sequential Matrix Multiplication
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
