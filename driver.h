#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define true  1
#define false 0

#define n_runs 1

double ** allocate_matrix(int dim);
int z_order_lookup(int x, int y);
struct indeces z_order_inverse_lookup(int value);
