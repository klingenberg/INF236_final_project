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
    
    double **C;
    double **A11, **A21, **A12, **A22;
    double **B11, **B21, **B12, **B22;
    double **C11, **C21, **C12, **C22;

    // help variables
    double **P1, **P2, **P3, **P4, **P5, **P6, **P7;
    double **S1, **S2, **S3, **S4, **S5, **S6, **S7, **S8;
    double **T1, **T2, **T3;
    
    // *********************************
    // Sequential Strassen Algorithm
    // *********************************
    
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

    // http://ftp.demec.ufpr.br/CFD/bibliografia/Higham_2002_Accuracy%20and%20Stability%20of%20Numerical%20Algorithms.pdf
    // page 436 (465 in pdf)
    // Winograd's form of Strassen algorithm, only 15 additions, not 18
    
    S1 = add(A21, A22, k);
    S2 = sub(S1, A11, k);
    S3 = sub(A11, A21, k);
    S4 = sub(A12, S2, k);
    S5 = sub(B12, B11, k);
    S6 = sub(B22, S5, k);
    S7 = sub(B22, B12, k);
    S8 = sub(S6, B21, k);
    
    P1 = sequential_strassen(S2, S6, k);
    P2 = sequential_strassen(A11, B11, k);
    P3 = sequential_strassen(A12, B21, k);
    P4 = sequential_strassen(S3, S7, k);
    P5 = sequential_strassen(S1, S5, k);
    P6 = sequential_strassen(S4, B22, k);
    P7 = sequential_strassen(A22, S8, k);
    
    T1 = add(P1, P2, k);
    T2 = add(T1, P4, k);
    T3 = add(T1, P5, k);
    
    C11 = add(P2, P3, k);
    C12 = add(T3, P6, k);
    C21 = sub(T2, P7, k);
    C22 = add(T2, P5, k);
    
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
    free(S1);
    free(S2);
    free(S3);
    free(S4);
    free(S5);
    free(S6);
    free(S7);
    free(S8);
    free(T1);
    free(T2);
    free(T3);

    return C;
}
