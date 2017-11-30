#include <stdio.h>
#include <stdlib.h>
#include "twb-quad.h"
#include "problem-spec.h"
#include "mesh.h"

double integrate_over_triangle(struct elem *ep, struct TWB_qdat *qdat,
		double (*f)(double x, double y))
{
	double x[3], y[3];
	double sum = 0.0;

	for (int i = 0; i < 3; i++) {
		x[i] = ep->n[i]->x;
		y[i] = ep->n[i]->y;
	}

	while (qdat->weight != -1.0) {
		double lambda1 = qdat->lambda1;
		double lambda2 = qdat->lambda2;
		double lambda3 = qdat->lambda3;
		double X = lambda1*x[0] + lambda2*x[1] + lambda3*x[2];
		double Y = lambda1*y[0] + lambda2*y[1] + lambda3*y[2];
		sum += qdat->weight * f(X, Y);
		qdat++;
	}

	sum *= ep->area / TWB_STANDARD_AREA;
	
	return sum;
}

void do_demo(struct problem_spec *spec, double a,
		struct TWB_qdat *qdat, char *filename)
{
	struct mesh *mesh = make_mesh(spec, a);
	printf("mesh of %d nodes, %d edges, %d elements\n",
			mesh->nnodes, mesh->nedges, mesh->nelems);

	double sum = 0.0;
	for (int i = 0; i < mesh->nelems; i++)
		sum +=integrate_over_triangle(&mesh->elems[i], qdat, spec->f);

	printf("the integral is %.12f\n", sum);

	free_mesh(mesh);
}

void show_usage_and_exit(char *progname)
{
	fprintf(stderr, "Usage: %s d a\n", progname);
	fprintf(stderr, "    d = twb integraton strength\n");
	fprintf(stderr, "    a = maximal triangle area\n");
	exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
	int d, n;
	char *endptr;
	double a;
	struct problem_spec *spec;
	struct problem_spec *triangle_with_hole(void);
	struct problem_spec *square(void);
	struct problem_spec *three_holes(int n);

	if (argc != 3)
		show_usage_and_exit(argv[0]);

	d = strtol(argv[1], &endptr, 10);
	if (*endptr != '\0')
		show_usage_and_exit(argv[0]);

	a = strtod(argv[2], &endptr);
	if (*endptr != '\0' || a <= 0.0)
		show_usage_and_exit(argv[0]);

	struct TWB_qdat *qdat = twb_qdat(&d, &n);

	printf("integrating with strength %d and %d points\n", d, n);

	spec = square();
	do_demo(spec, a, qdat, "square.gv");
	free_spec(spec);
	
	spec = three_holes(20);
	do_demo(spec, a, qdat, "three-holes.gv");
	free_spec(spec);

	return 0;
}
