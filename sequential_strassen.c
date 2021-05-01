void sub(double *C, double *A, double *B, int n) {
    int i;
    for(i = 0; i < n; i++) C[i] = A[i] - B[i];
}

void add(double *C, double *A, double *B, int n) {
    int i;
    for(i = 0; i < n; i++) C[i] = A[i] + B[i];
}

double * reorder_matrix(double ** X, int dim) {
    double *T;
    T = (double*) malloc(sizeof(double) * dim * dim);

    if (T == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }

    int i,j;
    for (i = 0; i < dim; i++) for (j = 0; j < dim; j++) T[z_order_lookup(j, i)] = X[i][j];

    return T;
}

void reorder_to_standard_matrix(double **C, double *X, int dim) {
    int i,j;
    for (i = 0; i < dim; i++) for (j = 0; j < dim; j++)  C[i][j] = X[z_order_lookup(j, i)];
}

double * allocate_1d_matrix(int dim) {
    double *T;
    T = (double*) malloc(sizeof(double) * dim * dim);

    if (T == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }

    return T;
}

void matmul(double * C,double * A, double * B, int dim) {
    
    int i, j, k;
    
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            C[z_order_lookup(j, i)] = 0.0;
            for(k = 0; k < dim; k++) {
                C[z_order_lookup(j, i)] += A[z_order_lookup(k, i)] * B[z_order_lookup(j, k)];
            } // k
        } // j
    } // i
}


int sequential_strassen_recursion(double *C, double *A, double *B, int n, double *X, double *Y){

    double *A11, *A21, *A12, *A22;
    double *B11, *B21, *B12, *B22;
    double *N1, *N2, *N3, *N4, *N5, *N6, *X_small, *Y_small;

    
    // *********************************
    // Sequential Strassen Algorithm
    // *********************************
    
    // Depth level less than  8:
    if (n <= 8) {
        matmul(C, A, B, n);
        return 0;
    }
    
    int k = n / 2;
    int kk = k*k;

    A11 = &A[0];
    A12 = &A[kk];
    A21 = &A[2*kk];
    A22 = &A[3*kk];

    B11 = &B[0];
    B12 = &B[kk];
    B21 = &B[2*kk];
    B22 = &B[3*kk];

    N1 = &X[0];
    N2 = &X[kk];
    N3 = &X[2*kk];
    X_small = &X[3*kk];

    N4 = &Y[0];
    N5 = &Y[kk];
    N6 = &Y[2*kk];
    Y_small = &Y[3*kk];
 

    // https://arxiv.org/pdf/0707.2347.pdf
    // Winograd's form of Strassen algorithm, only 15 additions, not 18
    // reduce number of temporary matrices to 6

    // S3 = A11 - A21
    sub(N1, A11, A21, kk);
    // T3 = B22 - B12
    sub(N2, B22, B12, kk);
    // P7 = S3 * T3
    sequential_strassen_recursion(N5, N1, N2, k, X_small, Y_small);
    
    // S1 = A21 + A22
    add(N1, A21, A22, kk);
    // T1 = B12 − B11
    sub(N2, B12, B11, kk);
    // P5 = S1 * T1
    sequential_strassen_recursion(N6, N1, N2, k, X_small, Y_small);

    // S2 = S1 − A11
    sub(N1, N1, A11, kk);
    // T2 = B22 − T1
    sub(N2, B22, N2, kk);
    // P6 = S2 * T2
    sequential_strassen_recursion(N4, N1, N2, k, X_small, Y_small);
    
    // S4 = A12 − S2
    sub(N1, A12, N1, kk);
    // P3 = S4 * B22
    sequential_strassen_recursion(N3, N1, B22, k, X_small, Y_small);

    // P1 = A11 * B11 
    sequential_strassen_recursion(N1, A11, B11, k, X_small, Y_small);

    // U2 = P1 + P6
    add(N4, N1, N4, kk);
    // U3 = U2 + P7
    add(N5, N4, N5, kk);
    // U4 = U2 + P5
    add(N4, N4, N6, kk);

    // U7 = U3 + P5
    add(N6, N5, N6, kk); // final C22
    // U5 = U4 + P3
    add(N4, N4, N3, kk); // final C12

    // T4 = T2 − B21
    sub(N2, N2, B21, kk);
    // P4 = A22 * T4
    sequential_strassen_recursion(N3, A22, N2, k, X_small, Y_small);

    // U6 = U3 − P4
    sub(N5, N5, N3, kk); // final C21

    // P2 = A12 * B21
    sequential_strassen_recursion(N3, A12, B21, k, X_small, Y_small);
    // U1 = P1 + P2
    add(N3, N1, N3, kk); // final C11

    int i,j;

    for(i = 0; i < kk; i++) {
        C[i] = N3[i];
        C[i + kk] = N4[i];
        C[i + 2*kk] = N5[i];
        C[i + 3*kk] = N6[i];
    }
/*
    printf("%d\n",n);

    printf("A:\n");
    for(i = 0; i < 4*kk; i++) {
        printf("%f ",A[i]);
    }
    printf("\n");
    for(i = 0; i < kk; i++) {
        printf("%f ",A11[i]);
    }
    printf("\n");
    for(i = 0; i < kk; i++) {
        printf("%f ",A12[i]);
    }
    printf("\n");
    for(i = 0; i < kk; i++) {
        printf("%f ",A21[i]);
    }
    printf("\n");
    for(i = 0; i < kk; i++) {
        printf("%f ",A22[i]);
    }
    printf("\n");
    printf("B:\n");
    for(i = 0; i < 4*kk; i++) {
        printf("%f ",B[i]);
    }
    printf("\n");
    for(i = 0; i < 4*kk; i++) {
        printf("%f ",C[i]);
    }
    printf("\n");
*/
    
    return 0;

}

double ** sequential_strassen(double **A, double **B, int n){
    double *X;
    X = allocate_1d_matrix(n);

    double *rA, *rB;
    rA = reorder_matrix(A, n);
    rB = reorder_matrix(B, n);

    // help variables
    double *N1, *N2;

    N1 = allocate_1d_matrix(n);
    N2 = allocate_1d_matrix(n);

    sequential_strassen_recursion(X, rA, rB, n, N1, N2);

    double **C;
    C = allocate_matrix(n);
    reorder_to_standard_matrix(C, X, n);

    free(X);
    free(rA);
    free(rB);
    free(N1);
    free(N2);

    return C;

}

