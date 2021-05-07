#include "morton_order.c"

void print_array_as_matrix(double *A,int dim) {
    int j;
    for (j = 0; j < dim*dim; j++) {
        if (j % dim == 0) printf("\n");
        printf("%.0f ", A[j]);
    }
    printf("\n");
}

void sub(double *C, double *A, double *B, int n) {
    int i;
    for(i = 0; i < n; i++) C[i] = A[i] - B[i];
}

void add(double *C, double *A, double *B, int n) {
    int i;
    for(i = 0; i < n; i++) C[i] = A[i] + B[i];
}

/*
double * reorder_to_z_array(double ** X, int dim) {
    double *T;
    T = allocate_array(dim*dim);

    int i,j;
    for (i = 0; i < dim; i++) for (j = 0; j < dim; j++) T[z_order_lookup(j, i)] = X[i][j];

    return T;
}

void reorder_to_standard_matrix(double **C, double *X, int dim) {
    int i,j;
    for (i = 0; i < dim; i++) for (j = 0; j < dim; j++)  C[i][j] = X[z_order_lookup(j, i)];
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


int * mdb_sequence(int n) {
    int i;
    int *S;
    S = (int*) malloc(sizeof(int) * n);

    S[0] = 0;

    for(i = 1; i < n; i++) {
        S[i] = (S[i-1] + 0xaaaaaaab) & 0x55555555;
    }

    return S;
}

int matmul_hardcoded_4x4(double * C,double * A, double * B, int n) {

    C[0] = A[0] * B[0] + A[1] * B[2] + A[4] * B[8] + A[5] * B[10];
    C[1] = A[0] * B[1] + A[1] * B[3] + A[4] * B[9] + A[5] * B[11];
    C[2] = A[2] * B[0] + A[3] * B[2] + A[6] * B[8] + A[7] * B[10];
    C[3] = A[2] * B[1] + A[3] * B[3] + A[6] * B[9] + A[7] * B[11];

    C[4] = A[0] * B[4] + A[1] * B[6] + A[4] * B[12] + A[5] * B[14];
    C[5] = A[0] * B[5] + A[1] * B[7] + A[4] * B[13] + A[5] * B[15];
    C[6] = A[2] * B[4] + A[3] * B[6] + A[6] * B[12] + A[7] * B[14];
    C[7] = A[2] * B[5] + A[3] * B[7] + A[6] * B[13] + A[7] * B[15];
    
    C[8] = A[8] * B[0] + A[9] * B[2] + A[12] * B[8] + A[13] * B[10];
    C[9] = A[8] * B[1] + A[9] * B[3] + A[12] * B[9] + A[13] * B[11];
    C[10] = A[10] * B[0] + A[11] * B[2] + A[14] * B[8] + A[15] * B[10];
    C[11] = A[10] * B[1] + A[11] * B[3] + A[14] * B[9] + A[15] * B[11];
    
    C[12] = A[8] * B[4] + A[9] * B[6] + A[12] * B[12] + A[13] * B[14];
    C[13] = A[8] * B[5] + A[9] * B[7] + A[12] * B[13] + A[13] * B[15];
    C[14] = A[10] * B[4] + A[11] * B[6] + A[14] * B[12] + A[15] * B[14];
    C[15] = A[10] * B[5] + A[11] * B[7] + A[14] * B[13] + A[15] * B[15];

    return 0;
}

int matmul_recursive(double * C,double * A, double * B, int n, double *H) {

    double *A11, *A21, *A12, *A22;
    double *B11, *B21, *B12, *B22;
    double *C11, *C21, *C12, *C22;
    double *H1, *H2;
    
    if (n == 1) {
        C[0] = A[0] * B[0];
        
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

    C11 = &C[0];
    C12 = &C[kk];
    C21 = &C[2*kk];
    C22 = &C[3*kk];

    H1 = &H[0];
    H2 = &H[kk];
    
    matmul_recursive(H1,A11,B11,k,H2);
    matmul_recursive(C11,A12,B21,k,H2);
    add(C11,H1,C11,kk);

    matmul_recursive(H1,A11,B12,k,H2);
    matmul_recursive(C12,A12,B22,k,H2);
    add(C12,H1,C12,kk);

    matmul_recursive(H1,A21,B11,k,H2);
    matmul_recursive(C21,A22,B21,k,H2);
    add(C21,H1,C21,kk);

    matmul_recursive(H1,A21,B12,k,H2);
    matmul_recursive(C22,A22,B22,k,H2);
    add(C22,H1,C22,kk);
    
    return 0;
}

void matmul_morton(double * C,double * A, double * B, int dim, int *S) {
    
    int i, j, k;

    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            C[2*S[i]+S[j]] = 0.0;
            for(k = 0; k < dim; k++) {
                C[2*S[i]+S[j]] += A[2*S[i]+S[k]] * B[2*S[k]+S[j]];
            } // k
        } // j
    } // i
}
*/

int matmul(double *C, double * A, double * B, int dim) {
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
}

int sequential_strassen_recursion(double *C, double *A, double *B, int n, double *X, int depth){

    double *A11, *A21, *A12, *A22;
    double *B11, *B21, *B12, *B22;
    double *N1, *N2, *N3, *N4, *N5, *N6, *X_small;

    
    // *********************************
    // Sequential Strassen Algorithm
    // *********************************
    
    // Depth level:
    if (n <= depth) { // depth) {
        //matmul_recursive(C, A, B, n, X);
        //matmul_morton(C, A, B, n, M);
        //matmul_hardcoded_4x4(C, A, B, n);
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
    X_small = &X[2*kk];

    N3 = &C[0];
    N4 = &C[kk];
    N5 = &C[2*kk];
    N6 = &C[3*kk];
 

    // https://arxiv.org/pdf/0707.2347.pdf
    // Winograd's form of Strassen algorithm, only 15 additions, not 18
    // 6 temporary matrices

    // S3 = A11 - A21
    sub(N1, A11, A21, kk);
    /*
    printf("S3 = A11 - A21\n");
    print_array_as_matrix(A11,k);
    print_array_as_matrix(A21,k);
    print_array_as_matrix(N1,k);
    */
    // T3 = B22 - B12
    sub(N2, B22, B12, kk);
    /*
    printf("T3 = B22 - B12\n");
    print_array_as_matrix(B22,k);
    print_array_as_matrix(B12,k);
    print_array_as_matrix(N2,k);
    */
    // P7 = S3 * T3
    sequential_strassen_recursion(N5, N1, N2, k, X_small, depth);
    /*
    printf("P7 = S3 * T3\n");
    print_array_as_matrix(N5,k);
    */
    
    // S1 = A21 + A22
    add(N1, A21, A22, kk);
    // T1 = B12 − B11
    sub(N2, B12, B11, kk);
    // P5 = S1 * T1
    sequential_strassen_recursion(N6, N1, N2, k, X_small, depth);

    // S2 = S1 − A11
    sub(N1, N1, A11, kk);
    // T2 = B22 − T1
    sub(N2, B22, N2, kk);
    // P6 = S2 * T2
    sequential_strassen_recursion(N4, N1, N2, k, X_small, depth);
    
    // S4 = A12 − S2
    sub(N1, A12, N1, kk);
    // P3 = S4 * B22
    sequential_strassen_recursion(N3, N1, B22, k, X_small, depth);

    // P1 = A11 * B11 
    sequential_strassen_recursion(N1, A11, B11, k, X_small, depth);

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
    sequential_strassen_recursion(N3, A22, N2, k, X_small, depth);

    // U6 = U3 − P4
    sub(N5, N5, N3, kk); // final C21

    // P2 = A12 * B21
    sequential_strassen_recursion(N3, A12, B21, k, X_small, depth);
    // U1 = P1 + P2
    add(N3, N1, N3, kk); // final C11

    return 0;
}

double ** parallel_strassen_recursion(double **A, double **B, int n, float *t){
    double mt1, mt2; // Timing variables

    double *R;
    R = allocate_array(n*n);

    int depth = 32;

    double *rA, *rB;
    rA = reorder_to_morton_array(A, n, depth);
    rB = reorder_to_morton_array(B, n, depth);

    //print_array_as_matrix(rA,n);
    //print_array_as_matrix(rB,n);

    // help variables
    double *H;

    H = allocate_array(3*(n*n)/4); // size 3/4 of original matrix

    mt1 = omp_get_wtime();
    sequential_strassen_recursion(R, rA, rB, n, H, depth);
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

