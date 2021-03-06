void parallel_sub(double *C, double *A, double *B, int n) {
    int i;
    
    #pragma omp for
    for(i = 0; i < n; i++) C[i] = A[i] - B[i];
}

void parallel_add(double *C, double *A, double *B, int n) {
    int i;
    
    #pragma omp for
    for(i = 0; i < n; i++) C[i] = A[i] + B[i];
}

void parallel_matmul_strassen(double *C, double * A, double * B, int dim) {
    int i, j, k;

    #pragma omp for private(j,k)
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            C[i*dim+j] = 0.0;
        } // j
        for(k = 0; k < dim; k++) {
            for (j = 0; j < dim; j++) {
                C[i*dim+j] += A[i*dim+k] * B[k*dim+j];
            } // j
        } // k
    } // i
}

int parallel_strassen_end(double *C, double *A, double *B, int n, double *X, int depth){

    double *A11, *A21, *A12, *A22;
    double *B11, *B21, *B12, *B22;
    double *N1, *N2, *N3, *N4, *N5, *N6, *X_small;
    
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
    X_small = &X[2*kk];

    N3 = &C[0];
    N4 = &C[kk];
    N5 = &C[2*kk];
    N6 = &C[3*kk];

    #pragma omp parallel
    {

    // S3 = A11 - A21
    parallel_sub(N1, A11, A21, kk);
    // T3 = B22 - B12
    parallel_sub(N2, B22, B12, kk);

    // P7 = S3 * T3
    parallel_matmul_strassen(N5, N1, N2, k);
    
    // S1 = A21 + A22
    parallel_add(N1, A21, A22, kk);
    
    // T1 = B12 − B11
    parallel_sub(N2, B12, B11, kk);

    // P5 = S1 * T1
    parallel_matmul_strassen(N6, N1, N2, k);

    // S2 = S1 − A11
    parallel_sub(N1, N1, A11, kk);
    
    // T2 = B22 − T1
    parallel_sub(N2, B22, N2, kk);
    
    // P6 = S2 * T2
    parallel_matmul_strassen(N4, N1, N2, k);
    
    // S4 = A12 − S2
    parallel_sub(N1, A12, N1, kk);
    
    // P3 = S4 * B22
    parallel_matmul_strassen(N3, N1, B22, k);

    // P1 = A11 * B11
    parallel_matmul_strassen(N1, A11, B11, k);

    // U2 = P1 + P6
    parallel_add(N4, N1, N4, kk);
    
    // U3 = U2 + P7
    parallel_add(N5, N4, N5, kk);

    // U4 = U2 + P5
    parallel_add(N4, N4, N6, kk);

    // U7 = U3 + P5
    parallel_add(N6, N5, N6, kk); // final C22
    
    // U5 = U4 + P3
    parallel_add(N4, N4, N3, kk); // final C12

    // T4 = T2 − B21
    parallel_sub(N2, N2, B21, kk);

    // P4 = A22 * T4
    parallel_matmul_strassen(N3, A22, N2, k);

    // U6 = U3 − P4
    parallel_sub(N5, N5, N3, kk); // final C21

    // P2 = A12 * B21
    parallel_matmul_strassen(N3, A12, B21, k);

    // U1 = P1 + P2
    parallel_add(N3, N1, N3, kk); // final C11
  
    }

    return 0;
}

int parallel_strassen_recursion(double *C, double *A, double *B, int n, double *X, int depth){

    double *A11, *A21, *A12, *A22;
    double *B11, *B21, *B12, *B22;
    double *N1, *N2, *N3, *N4, *N5, *N6, *X_small;
    
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
    X_small = &X[2*kk];

    N3 = &C[0];
    N4 = &C[kk];
    N5 = &C[2*kk];
    N6 = &C[3*kk];

    #pragma omp parallel
    {

    // S3 = A11 - A21
    parallel_sub(N1, A11, A21, kk);
    // T3 = B22 - B12
    parallel_sub(N2, B22, B12, kk);
    }
    // P7 = S3 * T3
    if (k == depth) {
        parallel_strassen_end(N5, N1, N2, k, X_small, depth);
    } else {
        parallel_strassen_recursion(N5, N1, N2, k, X_small, depth);
    }
    
    
    #pragma omp parallel
    {
    // S1 = A21 + A22
    parallel_add(N1, A21, A22, kk);
    
    // T1 = B12 − B11
    parallel_sub(N2, B12, B11, kk);
    }

    // P5 = S1 * T1
    if (k == depth) {
        parallel_strassen_end(N6, N1, N2, k, X_small, depth);
    } else {
        parallel_strassen_recursion(N6, N1, N2, k, X_small, depth);
    }

    #pragma omp parallel
    {
    // S2 = S1 − A11
    parallel_sub(N1, N1, A11, kk);
    
    // T2 = B22 − T1
    parallel_sub(N2, B22, N2, kk);
    }

    // P6 = S2 * T2
    if (k == depth) {
        parallel_strassen_end(N4, N1, N2, k, X_small, depth);
    } else {
        parallel_strassen_recursion(N4, N1, N2, k, X_small, depth);
    }
    
    
    #pragma omp parallel
    {
    // S4 = A12 − S2
    parallel_sub(N1, A12, N1, kk);
    }
    // P3 = S4 * B22
    if (k == depth) {
        parallel_strassen_end(N3, N1, B22, k, X_small, depth);
    } else {
        parallel_strassen_recursion(N3, N1, B22, k, X_small, depth);
    }

    // P1 = A11 * B11
    if (k == depth) {
        parallel_strassen_end(N1, A11, B11, k, X_small, depth);
    } else {
        parallel_strassen_recursion(N1, A11, B11, k, X_small, depth);
    }

    #pragma omp parallel
    {
    // U2 = P1 + P6
    parallel_add(N4, N1, N4, kk);
    
    // U3 = U2 + P7
    parallel_add(N5, N4, N5, kk);

    // U4 = U2 + P5
    parallel_add(N4, N4, N6, kk);

    // U7 = U3 + P5
    parallel_add(N6, N5, N6, kk); // final C22
    
    // U5 = U4 + P3
    parallel_add(N4, N4, N3, kk); // final C12

    // T4 = T2 − B21
    parallel_sub(N2, N2, B21, kk);
    }

    // P4 = A22 * T4
    if (k == depth) {
        parallel_strassen_end(N3, A22, N2, k, X_small, depth);
    } else {
        parallel_strassen_recursion(N3, A22, N2, k, X_small, depth);
    }

    #pragma omp parallel
    {
    // U6 = U3 − P4
    parallel_sub(N5, N5, N3, kk); // final C21
    }

    // P2 = A12 * B21
    if (k == depth) {
        parallel_strassen_end(N3, A12, B21, k, X_small, depth);
    } else {
        parallel_strassen_recursion(N3, A12, B21, k, X_small, depth);
    }

    #pragma omp parallel
    {
    // U1 = P1 + P2
    parallel_add(N3, N1, N3, kk); // final C11
    }

    return 0;
}

int parallel_strassen(double **C, double **A, double **B, int n, float *t){
    double mt1, mt2; // Timing variables

    //int layers = 3;
    //int depth = n/(1 << (layers-1));
    int depth = 512;
    int not_reordered_submatrix_size = depth/2;

    if (n<4*depth) {
        mt1 = omp_get_wtime();
        parallel_matmul(C, A, B, n);
        mt2 = omp_get_wtime();

        *t = mt2 - mt1;
        return 0;
    }

    double *R;
    R = allocate_array(n*n);

    double *rA, *rB;
    rA = reorder_to_morton_array(A, n, not_reordered_submatrix_size);
    rB = reorder_to_morton_array(B, n, not_reordered_submatrix_size);
    // help variables
    double *H;

    H = allocate_array(3*(n*n)/4);

    mt1 = omp_get_wtime();

    parallel_strassen_recursion(R, rA, rB, n, H, depth);
    mt2 = omp_get_wtime();

    *t = mt2 - mt1;
    
    reorder_back_morton_array(C, R, n, not_reordered_submatrix_size);
    
    free(R);
    free(rA);
    free(rB);
    free(H);

    return 0;

}

