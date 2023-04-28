#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "func.c"
#include "mpi.h"

extern double func  (double t, double x);
extern double fi    (double x);
extern double ksi   (double t);

extern double t_max;
extern double x_max;
extern double t_step;
extern double x_step;


void init_matrix(double **matrix) {


}

void fill_matrix(double **matrix, int rank, int size) {
	int i, j;
	int n_k = (int) (t_max / t_step) + 1;
        int n_m = (int) (x_max / x_step) + 1;
	double frac = t_step / x_step;
	int min_col = (n_m / size) * rank;
	int max_col = (n_m / size) * (rank + 1) - 1;
	MPI_Status status;

	
	if (rank == 0)
		min_col = 1;

	if (rank == size - 1)
		max_col = n_m - 1;

	printf("rank: %d size: %d min_col: %d max_col: %d\n", rank, size, min_col, max_col);

	if (rank == 0) {
		for (i = 0; i < n_k; i++) {
			matrix[i][0] = ksi(i * t_step);
		}
	}
/*
	for (j = min_col; j < max_col; j++) {
		matrix[0][j] = fi(j * x_step);
	}
*/

	for (i = 0; i < n_k; i++) {
		for (j = min_col; j <= max_col; j++) {			
			if (i == 0) {
				matrix[i][j] = fi(j * x_step);
			}

			if ((rank != 0) && (j == min_col)) {
                		MPI_Recv (&matrix[i][j - 1], 1, MPI_DOUBLE, rank - 1, 5, MPI_COMM_WORLD, &status);
        		}

			if (i < n_k - 1) {
				matrix[i + 1][j] = (2 * t_step * func((i + 0.5) * t_step, (j + 0.5) * x_step) / (1 + frac)) -
					(matrix[i + 1][j - 1] * (1 - frac) / (1 + frac)) + 
					(matrix[i][j - 1]) + (matrix[i][j] * (1 - frac) / (1 + frac));
			}
		}

		if (rank != (size - 1)) {
         	       MPI_Send (&matrix[i][max_col], 1, MPI_DOUBLE, rank + 1, 5, MPI_COMM_WORLD);
        	}
	}
}

int main( int argc, char **argv ){
	int i, j;
	int rank, size;
	int n_k = (int) (t_max / t_step) + 1;
        int n_m = (int) (x_max / x_step) + 1;

	MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

	double** matrix = (double **)calloc(n_k, sizeof(double *));

	for (i = 0; i < n_k; i++) {
		matrix[i] = (double *)calloc(n_m, sizeof(double));
	}

	printf("n_k: %d n_m: %d\n", n_k, n_m);

	fill_matrix(matrix, rank, size);

	for (i = 0; i < n_k; i++) {
		for(j = 0; j < n_m; j++) {
			printf("%10.2f", matrix[i][j]);
		}
		printf("\n");
	}

	for (i = 0; i < n_k; i++) {
		free(matrix[i]);
	}

	free(matrix);

	MPI_Finalize();

    	return 0;
}
