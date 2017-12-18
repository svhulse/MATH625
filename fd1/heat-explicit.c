#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"
#include "problem-spec.h"

static double get_error(struct problem_spec *spec,
	double *u, int n, double T)
{
	double err = 0.0;
	for (int j = 0; j < n+2; j++) {
		double x = spec->a + (spec->b - spec->a)/(n+1)*j;
		double diff = fabs(u[j] - spec->u_exact(x, T));
		if (diff > err)
			err = diff;
	}
	
	return err;
}

static void plot_curve(FILE *fp, double *u, int n, int steps, int k)
{
	for (int j = 0; j < n+2; j++)
		fprintf(fp, "%g %g %g\n",
			(double)k/steps, (double)j/(n+1), u[j]);
}

void show_usage(char *progname)
{
	fprintf(stderr, "Usage: %s d a\n", progname);
	fprintf(stderr, "    T : solve over 0 <= t <= T\n");
	fprintf(stderr, "    n : number of grid points\n");
	fprintf(stderr, "    s : number of time-slices\n");
	exit(EXIT_FAILURE);
}

static void heat_explicit(struct problem_spec *spec,
	double T, int n, int steps, char *gv_filename)
{
	FILE *fp;
	double *u, *v;
	double dx = (spec->b - spec->a)/(n+1);
	double dt = T/steps;
	double r = dt/(dx*dx);
	
	if (r > 0.5)
		printf("Warning: finite difference scheme unstable: r = %f\n", r);

	if ((fp = fopen(gv_filename, "w")) == NULL) {
		fprintf(stderr, "unable to open file '%s' for writing\n",
			gv_filename);
		return;
	}

	fprintf(fp, "# geomview script written by the function %s()\n",
		__func__);
	fprintf(fp, "{ appearance { +edge }\n");
	fprintf(fp, "MESH %d %d\n", n+2, steps+1);
	printf("%g < x < %g, 0 < t < %g, dx = %g, dt = %g, "
		"r = dt/dx^2 = %g\n",
		spec->a, spec->b, T, dx, dt, r);

	make_vector(u, n+2);
	make_vector(v, n+2);

	for (int j = 0; j < n+2; j++) {
		double x = spec->a + (spec->b - spec->a)/(n+1)*j;
		u[j] = spec->ic(x);
	}
	
	plot_curve(fp, u, n, steps, 0);

	for (int k = 1; k <= steps; k++) {
		double *tmp;

		for (int i = 0; i < n; i++)
			v[i+1] = r*u[i] + (1-2*r)*u[i+1] + r*u[i+2];

		double t = T*k/steps;
		v[0] = spec->bcL(t);
		v[n+1] = spec->bcR(t);

		tmp = v;
		v = u;
		u = tmp;

		plot_curve(fp, u, n, steps, k);
	}

	fprintf(fp, "}\n");
	fclose(fp);
	printf("geomview scripy written to file %s\n", gv_filename);
	if (spec->u_exact != NULL) {
		double err = get_error(spec, u, n, T);
		printf("max error at time %g is %g\n", T, err);
	}

	free_vector(u);
	free_vector(v);
}

int main(int argc, char **argv)
{
	struct problem_spec *heat1(void);
	struct problem_spec *heat2(void);
	struct problem_spec *heat3(void);
	struct problem_spec *heat4(void);
	char *endptr;
	double T;
	int n, steps;

	if (argc != 4)
		show_usage(argv[0]);

	T = strtod(argv[1], &endptr);
	if (*endptr != '\0' || T <= 0.0) 
		show_usage(argv[0]);

	n = strtol(argv[2], &endptr, 10);
	if (*endptr != '\0' || n < 1)
		show_usage(argv[0]);

	steps = strtol(argv[3], &endptr, 10);
	if (*endptr != '\0' || steps < 0)
		show_usage(argv[0]);

	heat_explicit(heat1(), T, n, steps, "ex1.gv");
	heat_explicit(heat2(), T, n, steps, "ex2.gv");
	heat_explicit(heat3(), T, n, steps, "ex3.gv");
	heat_explicit(heat4(), T, n, steps, "ex4.gv");

	return EXIT_SUCCESS;
}