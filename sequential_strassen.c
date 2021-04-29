void sub(int **C, int **A, int **B, int n) {
    int i, j;
    for(i = 0; i < n; i++) for(j = 0; j < n; j++) C[i][j] = A[i][j] - B[i][j];
}

void add(int **C, int **A, int **B, int n) {
    int i, j;
    for(i = 0; i < n; i++) for(j = 0; j < n; j++) C[i][j] = A[i][j] + B[i][j];
}

int ** allocate_submatrix_pointers(int i, int j, int n, int** M) {
    int** S;
    S = (int**) malloc(sizeof(int *) * n);

    int k;
    for(k = 0; k < n; k++) S[k] = M[k+i*n] + j * n;
    return S;
}

int sequential_strassen(int **C, int **A, int **B, int n){

    int **A11, **A21, **A12, **A22;
    int **B11, **B21, **B12, **B22;

    // help variables
    int **N1, **N2, **N3, **N4, **N5, **N6;

    
    // *********************************
    // Sequential Strassen Algorithm
    // *********************************
    
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return 0;
    }

    int k = n / 2;

    A11 = allocate_submatrix_pointers(0, 0, k, A);
    A21 = allocate_submatrix_pointers(1, 0, k, A);
    A12 = allocate_submatrix_pointers(0, 1, k, A);
    A22 = allocate_submatrix_pointers(1, 1, k, A);

    B11 = allocate_submatrix_pointers(0, 0, k, B);
    B21 = allocate_submatrix_pointers(1, 0, k, B);
    B12 = allocate_submatrix_pointers(0, 1, k, B);
    B22 = allocate_submatrix_pointers(1, 1, k, B);

    // https://arxiv.org/pdf/0707.2347.pdf
    // Winograd's form of Strassen algorithm, only 15 additions, not 18
    // reduce number of temporary matrices to 6

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

    int i,j;

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
