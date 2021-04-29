double sub(double **C, double **A, double **B, int n) {
    int i, j;
    for(i = 0; i < n; i++) for(j = 0; j < n; j++) C[i][j] = A[i][j] - B[i][j];
}

double add(double **C, double **A, double **B, int n) {
    int i, j;
    for(i = 0; i < n; i++) for(j = 0; j < n; j++) C[i][j] = A[i][j] + B[i][j];
}

int sequential_strassen(double **C, double **A, double **B, int n){

    double **A11, **A21, **A12, **A22;
    double **B11, **B21, **B12, **B22;

    // help variables
    double **N1, **N2, **N3, **N4, **N5, **N6;

    
    // *********************************
    // Sequential Strassen Algorithm
    // *********************************
    
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return 0;
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

    // reduce number of temporaries https://arxiv.org/pdf/0707.2347.pdf

    N1 = allocate_matrix(k);
    N2 = allocate_matrix(k);
    N3 = allocate_matrix(k); // C11
    N4 = allocate_matrix(k); // C12
    N5 = allocate_matrix(k); // C21
    N6 = allocate_matrix(k); // C22

    // S3 = A11 - A21
    sub(N1, A11, A21, k);
    // T3 = B22 - B12
    sub(N2, B22, B12, k);
    // P7 = S3 * T3
    sequential_strassen(N5, N1, N2, k);
    
    // S1 = A21 + A22
    add(N1, A21, A22, k);
    // T1 = B12 − B11
    sub(N2, B12, B11, k);
    // P5 = S1 * T1
    sequential_strassen(N6, N1, N2, k);

    // S2 = S1 − A11
    sub(N1, N1, A11, k);
    // T2 = B22 − T1
    sub(N2, B22, N2, k);
    // P6 = S2 * T2
    sequential_strassen(N4, N1, N2, k);
    
    // S4 = A12 − S2
    sub(N1, A12, N1, k);
    // P3 = S4 * B22
    sequential_strassen(N3, N1, B22, k);

    // P1 = A11 * B11 
    sequential_strassen(N1, A11, B11, k);

    // U2 = P1 + P6
    add(N4, N1, N4, k);
    // U3 = U2 + P7
    add(N5, N4, N5, k);
    // U4 = U2 + P5
    add(N4, N4, N6, k);

    // U7 = U3 + P5
    add(N6, N5, N6, k); // final C22
    // U5 = U4 + P3
    add(N4, N4, N3, k); // final C12

    // T4 = T2 − B21
    sub(N2, N2, B21, k);
    // P4 = A22 * T4
    sequential_strassen(N3, A22, N2, k);

    // U6 = U3 − P4
    sub(N5, N5, N3, k); // final C21

    // P2 = A12 * B21
    sequential_strassen(N3, A12, B21, k);
    // U1 = P1 + P2
    add(N3, N1, N3, k); // final C11

    for(i = 0; i < k; i++) {
        for(j = 0; j < k; j++) {
            C[i][j] = N3[i][j];
            C[i][j + k] = N4[i][j];
            C[k + i][j] = N5[i][j];
            C[k + i][k + j] = N6[i][j];
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
    
    free(N1);
    free(N2);
    free(N3);
    free(N4);
    free(N5);
    free(N6);

    return 0;

}
