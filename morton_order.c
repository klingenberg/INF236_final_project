// http://www.netlib.org/utk/people/JackDongarra/CCDSC-2016/slides/talk34-walker.pdf

double * allocate_array(int n) {
    double *T;
    T = (double*) malloc(sizeof(double) * n);

    if (T == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }

    return T;
}

void swap(double *xp,  double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void unshuffle(double *T, double *N, int n) {
    int i,j,k;
    int col,row;

    for (i=0;i<n/4;i++) { 
        for (j=0;j<n;j++) {
            for (k=0;k<2;k++) {
                
                row = 2*i; if (j>=n/2) row++;
                col = j % (n/2) + k * (n/2);

                N[row*n+col] = T[i*n+k*(n*n/4)+j];
            } // k
        } // j
    } // i

    for (i = 0; i < n*n/2; i++) { 
        T[i] = N[i];
    }
}

void shuffle(double *T, double *N, int n) {
    int i,j,k;
    int col,row;

    for (i=0;i<n/4;i++) { 
        for (j=0;j<n;j++) {
            for (k=0;k<2;k++) {
                
                row = 2*i; if (j>=n/2) row++;
                col = j % (n/2) + k * (n/2);

                N[i*n+k*(n*n/4)+j] = T[row*n+col];
            } // k
        } // j
    } // i

    for (i = 0; i < n*n/2; i++) { 
        T[i] = N[i];
    }
}

void mortonOrder(double *A, double *H, int n, int b) { 
    int m,submatrix;
    for(m = b; m <n;m = m *2) {
        int size = (m*m);
        for(submatrix = 0; submatrix< n*n;submatrix = submatrix+size*2) {
            unshuffle(A+submatrix,H,m*2);
        }
    }
}

void mortonOrderBack(double *A, double *H, int n, int b) {            
    int m,submatrix;
    for(m = n/2; m >=b;m = m /2) {
        int size = (m*m);
        for(submatrix = 0; submatrix< n*n; submatrix = submatrix+size*2) {
            shuffle(A+submatrix,H,m*2);
        }
    }
}

double * reorder_to_morton_array(double ** X, int n, int depth) {
    double *T;
    int dim = n*n;
    int i,j;
    T = allocate_array(dim);

    for(i=0;i<n;i++) for(j=0;j<n;j++) T[i*n+j] = X[i][j];

    double *help;
    help = allocate_array(dim);

    mortonOrderBack(T,help,n,depth);
    free(help);

    return T;
}

void reorder_back_morton_array(double ** X, double *T, int n, int depth) {

    int dim = n*n;
    int i,j;

    double *help;
    help = allocate_array(dim);

    mortonOrder(T,help,n,depth);

    free(help);

    for(i=0;i<n;i++) for(j=0;j<n;j++) X[i][j] = T[i*n+j];
}

/*
int main(int argc, char *argv[]) {
    
    double **M;
    int i, j;
    
    int dim = 16; 
    int d = 4;

    M = allocate_matrix(dim);

    printf("Matrix: \n");
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            M[i][j] = i*dim+j;
            printf("%.0f ",M[i][j]);
        }printf("\n");
    }

    double *R;
    R = reorder_to_z_array(M, dim, d);


    printf("\nMorton:");
    for (j = 0; j < dim*dim; j++) {
        if (j % dim == 0) printf("\n");
        printf("%.0f ", R[j]);
    }printf("\n");

    reorder_back_z_array(M, R, dim, d);

    printf("\n");
    printf("Reordered back:\n");
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            M[i][j] = i*dim+j;
            printf("%.0f ",M[i][j]);
        }printf("\n");
    }

}
*/