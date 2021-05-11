#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define true  1
#define false 0

#define n_runs 3

double ** allocate_matrix(int dim, double *ptr);
