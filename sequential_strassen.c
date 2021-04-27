double ** allocate_matrix(int dim) {
    double **C;
    double * ptr;
    C = malloc(sizeof(double *) * dim);
    ptr = malloc(sizeof(double) * dim * dim);
    
    int i;
    for(i = 0; i < dim; i++) C[i] = ptr + (i * dim);
    
    if (C == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }
    return C;
}

double ** sub(double **A, double **B, int n) {
    double ** C = allocate_matrix(n);
    int i, j;
    for(i = 0; i < n; i++) for(j = 0; j < n; j++) C[i][j] = A[i][j] - B[i][j];
    return C;
}

double ** add(double **A, double **B, int n) {
    double ** C = allocate_matrix(n);
    int i, j;
    for(i = 0; i < n; i++) for(j = 0; j < n; j++) C[i][j] = A[i][j] + B[i][j];
    return C;
}

double ** sequential_strassen(double **A, double **B, int n){
    // handling the base case
    if (n == 1) {
        double **C = allocate_matrix(1);
        C[0][0] = A[0][0] * B[0][0];
        return C;
    }
    int k = n / 2;
    
    double **A11, **A21, **A12, **A22, **B11, **B21, **B12, **B22, **C11, **C21, **C12, **C22, **P1, **P2, **P3, **P4, **P5, **P6, **P7;
    
    A11 = allocate_matrix(k);
    A21 = allocate_matrix(k);
    A12 = allocate_matrix(k);
    A22 = allocate_matrix(k);
    B11 = allocate_matrix(k);
    B21 = allocate_matrix(k);
    B12 = allocate_matrix(k);
    B22 = allocate_matrix(k);
    
    int i, j;
    
    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) {
            A11[i][j] = A[i][j];
            B11[i][j] = B[i][j];
            A21[i][j] = A[k + i][j];
            B21[i][j] = B[k + i][j];
            A12[i][j] = A[i][k + j];
            B12[i][j] = B[i][k + j];
            A22[i][j] = A[k + i][k + j];
            B22[i][j] = B[k + i][k + j];
        }
    }
    
    P1 = sequential_strassen(add(A11, A22, k), add(B11, B22, k), k);
    P2 = sequential_strassen(add(A21, A22, k), B11, k);
    P3 = sequential_strassen(A11, sub(B12, B22, k), k);
    P4 = sequential_strassen(A22, sub(B21, B11, k), k);
    P5 = sequential_strassen(add(A11, A12, k), B22, k);
    P6 = sequential_strassen(sub(A21, A11, k), add(B11, B12, k), k);
    P7 = sequential_strassen(sub(A12, A22, k), add(B21, B22, k), k);
    
    
    C11 = sub(add(add(P1, P4, k), P7, k), P5, k);
    C21 = add(P2, P4, k);
    C12 = add(P3, P5, k);
    C22 = sub(add(add(P1, P3, k), P6, k), P2, k);
    
    double ** C = allocate_matrix(n);

    for(i = 0; i < k; i++) {
        for(j = 0; j < k; j++) {
            C[i][j] = C11[i][j];
            C[k + i][j] = C21[i][j];
            C[i][j + k] = C12[i][j];
            C[k + i][k + j] = C22[i][j];
        }
    }
    
    free(A11);
    free(A21);
    free(A12);
    free(A22);
    free(B11);
    free(B21);
    free(B12);
    free(B22);
    free(C11);
    free(C21);
    free(C12);
    free(C22);
    free(P1);
    free(P2);
    free(P3);
    free(P4);
    free(P5);
    free(P6);
    free(P7);
    
    return C;
}
