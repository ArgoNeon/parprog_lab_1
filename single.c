#include <stdio.h>
#include "stdlib.h"
#include "func.c"
#include "mpi.h"

extern double func  (double t, double x);
extern double fi    (double x);
extern double ksi   (double t);

extern double t_max;
extern double x_max;
extern double t_step;
extern double x_step;

void write_to_csv(double** matrix, int n_i, int n_j) {
        int i, j;
        FILE *out = fopen("single_output.csv", "w+");

        for (i = 0; i < n_i; i++) {
                for (j = 0; j < n_j; j++) {
                        fprintf (out, "%10.2f", matrix[i][j]);
                if (j != n_j - 1)
                        fprintf (out, ",");
        }

                fprintf (out, "\n");
        }

        fclose(out);
}

void fill_matrix(double **matrix, int n_k, int n_m) {
	int i, j;
	double frac = t_step / x_step;

	for (j = 0; j < n_m; j++) {
		matrix[0][j] = fi(j * x_step);
	}

	for (i = 0; i < n_k; i++) {
		matrix[i][0] = ksi(i * t_step);
	}

	for (i = 0; i < n_k - 1; i++) {
		for (j = 0; j < n_m; j++) {
			matrix[i + 1][j] = (2 * t_step * func((i + 0.5) * t_step, (j + 0.5) * x_step) / (1 + frac)) -
				(matrix[i + 1][j - 1] * (1 - frac) / (1 + frac)) + 
				(matrix[i][j - 1]) + (matrix[i][j] * (1 - frac) / (1 + frac));
		}
	}
}

int main( int argc, char **argv ){
	int i;
	int n_k = (int) (t_max / t_step) + 1;
	int n_m = (int) (x_max / x_step) + 1;
	double** matrix = (double **)calloc(n_k, sizeof(double *));

	for (i = 0; i < n_k; i++) {
		matrix[i] = (double *)calloc(n_m, sizeof(double));
	}

	printf("n_k: %d n_m: %d\n", n_k, n_m);

	double start = MPI_Wtime();

	fill_matrix(matrix, n_k, n_m);

	double end = MPI_Wtime();

	printf ("time: %f seconds\n", end - start);

	write_to_csv(matrix, n_k, n_m);

/*	double end = MPI_Wtime();

        printf ("time: %f seconds\n", end - start);
*/

	for (i = 0; i < n_k; i++) {
                free(matrix[i]);
        }

	free(matrix);

    	return 0;
}
