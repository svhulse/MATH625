#include <stdio.h>
#include <suitesparse/umfpack.h>
#include "array.h"

#define check_status(status) if (status != UMFPACK_OK) \
					   exit_on_error(status, __FILE__, __LINE__)

void exit_on_error(int status, const char *file, int line)
{
	if (status == UMFPACK_ERROR_out_of_memory)
		fprintf(stderr, "*** file %s, line %d: "
				"UMFPACK out of memory\n", file, line);
	else if (status == UMFPACK_WARNING_singular_matrix)
		fprintf(stderr, "*** file %s, line %d: "
				"UMFPACK recieved singular matrix\n", file, line);
	else
		fprintf(stderr, "*** file %s, line %d: "
				"trouble in UMFPACK call\n", file, line);
}

int main(void)
{
	double *b;
	double *x;
	int n = 5;
	int m = 5;
	int N = n*m;

	double *Tx;
	int *Ti, *Tj;

	double *Ax;
	int *Ai, *Ap;

	void *Numeric;
	void *Symbolic;

	int status;

	make_vector(Tx, n * m);
	make_vector(Ti, n * m);
	make_vector(Tj, n * m);

	make_vector(Ax, N);
	make_vector(Ai, N);
	make_vector(Ap, n+1);
	
	make_vector(x, n);
	make_vector(b, n);
	b[0] = 1;
	b[1] = 2;
	b[2] = 3;
	b[3] = 4;
	b[4] = 5;

	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			Ti[(j*n)+i] = i;
			Tj[(j*n)+i] = j;
			Tx[(j*n)+i] = 1.0 / (i + j + 1);
		}
	}

	status = umfpack_di_triplet_to_col(n, n, N, Ti, Tj, Tx, Ap, Ai, Ax, NULL);
	check_status(status);
	status = umfpack_di_symbolic(n, m, Ap, Ai, Ax, &Symbolic, NULL, NULL);
	check_status(status);
	status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
	check_status(status);
	status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
	check_status(status);

	print_vector("%10.1f", x, n);

	umfpack_di_free_symbolic(&Symbolic);
	umfpack_di_free_numeric(&Numeric);

	free_vector(x);
	free_vector(b);
	free_vector(Ax);
	free_vector(Ai);
	free_vector(Ap);
	free_vector(Tx);
	free_vector(Ti);
	free_vector(Tj);

	return 0;	
}
