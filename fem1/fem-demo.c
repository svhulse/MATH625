#include <stdio.h>
#include <stdlib.h>
#include "poisson.h"
#include "mesh.h"
#include "twb-quad.h"
#include "plot-with-geomview.h"

static void do_demo(struct problem_spec *spec,
    double a, int d, char *gv_filename)
{
	struct mesh *mesh;
	mesh = make_mesh(spec, a);
	
	printf("nodes = %d, edges = %d, elems = %d\n",
		mesh->nnodes, mesh->nedges, mesh->nelems);
	
	poisson_solve(spec, mesh, d);

	if (spec->u_exact != NULL) {
		struct errors errs = eval_errors(spec, mesh, d);
		printf("errors: L^infty = %g, L^2 = %g, energy norm = %g\n",
		errs.Linfty, errs.L2norm, errs.energy);
    	}

    	plot_with_geomview_zhue(mesh, gv_filename);
    	free_mesh(mesh);
}

void show_usage(char *progname)
{
	fprintf(stderr, "Usage: %s d a\n", progname);
	fprintf(stderr, "    d = twb integraton strength\n");
	fprintf(stderr, "    a = maximal triangle area\n");
	exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
	int d;
	char *endptr;
	double a;
	struct problem_spec *spec;
	
	void free_spec(struct problem_spec *spec);
	struct problem_spec *triangle_with_hole(void);
	struct problem_spec *square(void);
	struct problem_spec *three_holes(int n);

	if (argc != 3)
		show_usage(argv[0]);

	d = strtol(argv[1], &endptr, 10);
	if (*endptr != '\0')
		show_usage(argv[0]);

	a = strtod(argv[2], &endptr);
	if (*endptr != '\0' || a <= 0.0)
		show_usage(argv[0]);
	
	spec = square();
	do_demo(spec, a, d, "triangle_with_hole.gv");
	//free_spec(spec);

	return 0;
}
