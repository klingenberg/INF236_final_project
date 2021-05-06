void parallel_sub(double *C, double *A, double *B, int n) {
    int i;
    
    for(i = 0; i < n; i++) C[i] = A[i] - B[i];
}

void parallel_add(double *C, double *A, double *B, int n) {
    int i;

    for(i = 0; i < n; i++) C[i] = A[i] + B[i];
}

double * parallel_matmul_strassen(double *C, double * A, double * B, int dim) {
    int i, j, k;
    
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

    return C;
}

int parallel_strassen_recursion(double *C, double *A, double *B, int n, double *X, int depth){

    double *A11, *A21, *A12, *A22;
    double *B11, *B21, *B12, *B22;
    double *N1, *N2, *N3, *N4, *N5, *N6, *X_small;

    
    // *********************************
    // Sequential Strassen Algorithm
    // *********************************
    
    // Depth level:
    if (n <= depth) {
        parallel_matmul_strassen(C, A, B, n);
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
    X_small = &X[2*kk];

    N3 = &C[0];
    N4 = &C[kk];
    N5 = &C[2*kk];
    N6 = &C[3*kk];
 

    // https://arxiv.org/pdf/0707.2347.pdf
    // Winograd's form of Strassen algorithm, only 15 additions, not 18
    // 6 temporary matrices

    
    
    // S3 = A11 - A21
//    #pragma omp task shared(N1, A11, A21, kk)
    parallel_sub(N1, A11, A21, kk);
    
    // T3 = B22 - B12
//    #pragma omp task shared(N2, B22, B12, kk)
    parallel_sub(N2, B22, B12, kk);
    
//    #pragma omp taskwait
    
    
    // P7 = S3 * T3
//    #pragma omp task shared(N5, N1, N2, k, X_small, depth)
    parallel_strassen_recursion(N5, N1, N2, k, X_small, depth);
    
//    #pragma omp taskwait
    
    // S1 = A21 + A22
    
//    #pragma omp task shared(N1, A21, A22, kk)
    parallel_add(N1, A21, A22, kk);
    
    // T1 = B12 − B11
//    #pragma omp task shared(N2, B12, B11, kk)
    parallel_sub(N2, B12, B11, kk);
    
//    #pragma omp taskwait
    
    // P5 = S1 * T1
//    #pragma omp task shared(N6, N1, N2, k, X_small, depth)
    parallel_strassen_recursion(N6, N1, N2, k, X_small, depth);
    
//    #pragma omp taskwait

    // S2 = S1 − A11
    
//    #pragma omp task shared(N1, A11, kk)
    parallel_sub(N1, N1, A11, kk);
    
    // T2 = B22 − T1
//    #pragma omp task shared(N2, B22, N2, kk)
    parallel_sub(N2, B22, N2, kk);
    
//    #pragma omp taskwait
    
    // P6 = S2 * T2
//    #pragma omp task shared(N4, N1, N2, k, X_small, depth)
    parallel_strassen_recursion(N4, N1, N2, k, X_small, depth);
    
//    #pragma omp taskwait
    
    // S4 = A12 − S2
//    #pragma omp task shared(N1, A12, kk)
    parallel_sub(N1, A12, N1, kk);
    
//    #pragma omp taskwait
    
    // P3 = S4 * B22
//    #pragma omp task shared(N3, N1, B22, k, X_small, depth)
    parallel_strassen_recursion(N3, N1, B22, k, X_small, depth);

    // P1 = A11 * B11
//    #pragma omp task shared(N1, A11, B11, k, X_small, depth)
    parallel_strassen_recursion(N1, A11, B11, k, X_small, depth);
    
//    #pragma omp taskwait

    // U2 = P1 + P6
//    #pragma omp task shared(N4, N1, kk)
    parallel_add(N4, N1, N4, kk);
    
//    #pragma omp taskwait
    
    // U3 = U2 + P7
//    #pragma omp task shared(N5, N4, kk)
    parallel_add(N5, N4, N5, kk);
    
//    #pragma omp taskwait
    
    // U4 = U2 + P5
//    #pragma omp task shared(N4, N6, kk)
    parallel_add(N4, N4, N6, kk);
    
//    #pragma omp taskwait

    // U7 = U3 + P5
//    #pragma omp task shared(N5, N6, kk)
    parallel_add(N6, N5, N6, kk); // final C22
    
    // U5 = U4 + P3
//    #pragma omp task shared(N4, N3, kk)
    parallel_add(N4, N4, N3, kk); // final C12
    
    // T4 = T2 − B21
//    #pragma omp task shared(N2, B21, kk)
    parallel_sub(N2, N2, B21, kk);
    
//    #pragma omp taskwait
    
    // P4 = A22 * T4
//    #pragma omp task shared(N3, A22, N2, k, X_small, depth)
    parallel_strassen_recursion(N3, A22, N2, k, X_small, depth);

//    #pragma omp taskwait

    // U6 = U3 − P4
//    #pragma omp task shared(N5, N3, kk)
    parallel_sub(N5, N5, N3, kk); // final C21
    
    // P2 = A12 * B21
//    #pragma omp task shared(N3, A12, B21, k, X_small, depth)
    parallel_strassen_recursion(N3, A12, B21, k, X_small, depth);

//    #pragma omp taskwait
    
    // U1 = P1 + P2
//    #pragma omp task shared(N1, N3, kk)
    parallel_add(N3, N1, N3, kk); // final C11

//    #pragma omp taskwait
    
    return 0;
}

double ** parallel_strassen(double **A, double **B, int n, float *t){
    double mt1, mt2; // Timing variables

    double *R;
    R = allocate_array(n*n);

    int depth = 32;

    double *rA, *rB;
    rA = reorder_to_morton_array(A, n, depth);
    rB = reorder_to_morton_array(B, n, depth);
    // help variables
    double *H;

    H = allocate_array(3*(n*n)/4); // size 3/4 of original matrix

    mt1 = omp_get_wtime();
    
    #pragma omp parallel
    #pragma omp single nowait
    parallel_strassen_recursion(R, rA, rB, n, H, depth);
    mt2 = omp_get_wtime();

    *t = mt2 - mt1;
    
    double **C = allocate_matrix(n);
    
    reorder_back_morton_array(C, R, n, depth);
    
    free(R);
    free(rA);
    free(rB);
    free(H);

    return C;

}

