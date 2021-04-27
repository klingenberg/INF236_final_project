double ** sequential_strassen(double ** A, double ** B, int dim) {
    
    int i, j, k;
    double **C;
    double * ptr;
    C = malloc(sizeof(double *)*dim);
    ptr = malloc(sizeof(double)*dim*dim);
    
    for(i = 0; i < dim; i++) {
        C[i] = ptr + (i*dim);
    }
    
    if (C == NULL) {
        printf("Unable to allocate memory, exiting \n");
        exit(0);
    }
    
    // *********************************
    // Sequential Strassen Algorithm
    // *********************************
    
    //TODO
    
    return C;
}

double ** sequential_strassen_rec(){}
