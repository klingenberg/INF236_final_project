// we used this as inspiration http://www.netlib.org/utk/people/JackDongarra/CCDSC-2016/slides/talk34-walker.pdf

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

    for (i=0;i<n/4;i++) {                           // go through one quarter of array length
        for (j=0;j<n;j++) {                         // go through array length
            for (k=0;k<2;k++) {                     // two iterations
                                                    // ^ that all together goes through the whole array
                
                row = 2*i; if (j>=n/2) row++;       // row position in new matrix
                col = j % (n/2) + k * (n/2);        // column position in new matrix

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

    for (i=0;i<n/4;i++) {                           // go through one quarter of array length
        for (j=0;j<n;j++) {                         // go through array length
            for (k=0;k<2;k++) {                     // two iterations
                                                    // ^ that all together goes through the whole array

                row = 2*i; if (j>=n/2) row++;       // row position in old matrix
                col = j % (n/2) + k * (n/2);        // column position in old matrix
                
                N[i*n+k*(n*n/4)+j] = T[row*n+col];
            } // k
        } // j
    } // i

    for (i = 0; i < n*n/2; i++) { 
        T[i] = N[i];
    }
}

void mortonOrderBack(double *A, double *H, int n, int b) { 
    int m,submatrix;
    for(m = b; m <n;m = m *2) {      // go through matrix sizes powers of 2
        int size = (m*m);
        for(submatrix = 0; submatrix< n*n;submatrix = submatrix+size*2) {
            unshuffle(A+submatrix,H,m*2);
        }
    }
}

void mortonOrder(double *A, double *H, int n, int b) {            
    int m,submatrix;
    for(m = n/2; m >=b;m = m /2) {     // go through matrix sizes powers of 2
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

    // change from 2D to 1D matrix
    for(i=0;i<n;i++) for(j=0;j<n;j++) T[i*n+j] = X[i][j];

    double *help;
    help = allocate_array(dim);

    mortonOrder(T,help,n,depth);
    free(help);

    return T;
}

void reorder_back_morton_array(double ** X, double *T, int n, int depth) {

    int dim = n*n;
    int i,j;

    double *help;
    help = allocate_array(dim);

    mortonOrderBack(T,help,n,depth);

    free(help);

    // change from 1D to 2D matrix
    for(i=0;i<n;i++) for(j=0;j<n;j++) X[i][j] = T[i*n+j];
}