#include <stdio.h>
#include <stdlib.h>
#include "poisson.h"
#include "mesh.h"

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