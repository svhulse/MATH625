#include <suitesparse/umfpack.h>
#include <stdio.h>
#include "array.h"

void exit_on_error(int status, const char *file, int line)
{
	if (status == UMFPACK_ERROR_out_of_memory)
		fprinf(stderr, "*** file %s, line %d: " 
				"UMFPACK out of memory\n", file, line);
	else if (status == UMFPACK_WARNING_singular_matrix)
		fprintf(stderr, "*** file %s, line, %d: "
				"UMFPACK recieved singular matrix\n", file, line);
	else
		fprintf(stderr, "*** file %s, line %d: "
				"trouble in UMFPACK call\n", file, line);
}

int main(void)
{
	int n = 5;
	int N = 20;
	double *Tx;
	int *Ti, *Tj;
	double *Ax;
	int *Ap, *Ai;
	int nz;
	void *Numeric;
	void *Symbolic;
	double *b, *x;
	int status;

	make_vector(Tx, N);
	make_vector(Ti, N);
	make_vector(Tj, N);

	Ti[0] = 0;
	Ti[1] = 1;
	Ti[2] = 0;
	Ti[3] = 
	Ti[4] = 
	Ti[5] =
	Ti[6] =
	Ti[7] = 
	Ti[8] = 
	Ti[9] = 
	Ti[10] =
	Ti[11] =

	Tj[0] = 
	Tj[1] = 
	Tj[2] = 
	Tj[3] = 
	Tj[4] = 
	Tj[5] = 
	Tj[6] = 
	Tj[7] = 
	Tj[8] = 
	Tj[9] = 
	Tj[10] = 
	Tj[11] = 

	Tx[0] =
	Tx[1] = 
	Tx[2] =
	Tx[3] = -1.0
	Tx[4] = 3.0
	Tx[5] = 4.0;
	Tx[6] = -3.0;
	Tx[7] = 1.0;
	Tx[8] = 2.0;
	Tx[9] = 2.0;
	Tx[10] = 6.0;
	Tx[11] = 1.0;

	make_vector(Ax, N);
	make_vector(Ai, N);
	make_vector(Ap, n+1); 

	status = mfpack_di_triplet_to_col(n, n, N, Ti, Tj, Tx, Ap, Ai, Ax, NULL);
	if (status != UMFPACK_OK)
		exit_on_error(status, __FILE__, __LINE__);
	
	nz = Ap[n];
	printf("matrix ");

	umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
	umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);

	make_vector(b, n);
	make_vector(x, n);
	b[0] = 8.0;
	b[1] = 45.0;
	b[2] =	-3.0;
	b[3] = 3.0;
	b[4] = 19.0;

	umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
	printf("b = ");
	print_vector("%g ", b, n);


	printf("x = ");
	print_vector("%g ", x, n);

	free_vector(Tx);
	free_vector(Ti);
	free_vector(Tj);
	free_vector(Ax);
	free_vector(Ai);
	free_vector(Ap);
	free_vector(b);
	free_vector(x);

	umfpack_di_free_symbolic(&Symbolic);
	umfpack_di_free_numeric(&Numeric);


	return 0;
}
