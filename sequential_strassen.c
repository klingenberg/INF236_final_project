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
    
    double **C, **A11, **A21, **A12, **A22, **B11, **B21, **B12, **B22, **C11, **C21, **C12, **C22, **P1, **P2, **P3, **P4, **P5, **P6, **P7, **A11pA22, **B11pB22, **A21pA22, **B12mB22, **B21mB11, **A11mA12, **A21mA11, **B11pB12, **A12pA22, **B21pB22, **P1pP4, **P1pP4pP7, **P1pP3, **P1pP3pP6;
    
    // handling the base case
    if (n == 1) {
        C = allocate_matrix(1);
        C[0][0] = A[0][0] * B[0][0];
        
        return C;
    }
    int k = n / 2;

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
    
    A11pA22 = add(A11, A22, k);
    B11pB22 = add(B11, B22, k);
    A21pA22 = add(A21, A22, k);
    B12mB22 = sub(B12, B22, k);
    B21mB11 = sub(B21, B11, k);
    A11mA12 = add(A11, A12, k);
    A21mA11 = sub(A21, A11, k);
    B11pB12 = add(B11, B12, k);
    A12pA22 = sub(A12, A22, k);
    B21pB22 = add(B21, B22, k);
    
    P1 = sequential_strassen(A11pA22, B11pB22, k);
    P2 = sequential_strassen(A21pA22, B11, k);
    P3 = sequential_strassen(A11, B12mB22, k);
    P4 = sequential_strassen(A22, B21mB11, k);
    P5 = sequential_strassen(A11mA12, B22, k);
    P6 = sequential_strassen(A21mA11, B11pB12, k);
    P7 = sequential_strassen(A12pA22, B21pB22, k);
    
    P1pP4 = add(P1, P4, k);
    P1pP4pP7 = add(P1pP4, P7, k);
    P1pP3 = add(P1, P3, k);
    P1pP3pP6 = add(P1pP3, P6, k);
    
    C11 = sub(P1pP4pP7, P5, k);
    C21 = add(P2, P4, k);
    C12 = add(P3, P5, k);
    C22 = sub(P1pP3pP6, P2, k);
    
    C = allocate_matrix(n);

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
    free(A11pA22);
    free(B11pB22);
    free(A21pA22);
    free(B12mB22);
    free(B21mB11);
    free(A11mA12);
    free(A21mA11);
    free(B11pB12);
    free(A12pA22);
    free(B21pB22);
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
    free(P1pP4);
    free(P1pP3);
    free(P1pP4pP7);
    free(P1pP3pP6);

    return C;
}
