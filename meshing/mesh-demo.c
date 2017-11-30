#include <stdio.h>
#include <stdlib.h>
#include "problem-spec.h"
#include "mesh-to-eps.h"

void show_usage(char *progname)
{
	fprintf(stderr, "Usage: %s a\n", progname);
	fprintf(stderr, "   a = maximal triangle area\n");
}

void generate_mesh(struct problem_spec *(*mesh_spec)(void), double a, char *name)
{
	struct problem_spec *spec;
	struct mesh *mesh;

	spec = (*mesh_spec)();
	mesh = make_mesh(spec, a);	
	printf("mesh with %d nodes, %d edges, %d elements\n",
			mesh->nnodes, mesh->nedges, mesh->nelems);
	mesh_to_eps(mesh, name);
	free_mesh(mesh);
}

void generate_mesh_1param(struct problem_spec *(*mesh_spec)(int), int n, double a, char *name)
{
	void free_spec(struct problem_spec *spec);
	struct problem_spec *spec;
	struct mesh *mesh;

	spec = (*mesh_spec)(n);
	mesh = make_mesh(spec, a);
	printf("mesh with %d nodes, %d edges, %d elements\n",
			mesh->nnodes, mesh->nedges, mesh->nelems);
	mesh_to_eps(mesh, name);

	free_spec(spec);
	free_mesh(mesh);
}

int main(int argc, char **argv)
{
	struct problem_spec *triangle_with_hole(void);
	struct problem_spec *square(void);
	struct problem_spec *three_holes(int n);
	struct problem_spec *annulus(int n);
	
	double a;
	char *endptr;

	if (argc != 2) {
		show_usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	a = strtod(argv[1], &endptr);
	if (*endptr != '\0' || a <= 0.0) {
		show_usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	generate_mesh(triangle_with_hole, a, "triangle-with-hole.eps");
	generate_mesh(square, a, "square.eps");
	generate_mesh_1param(annulus, 20, a, "annulus.eps");
	generate_mesh_1param(three_holes, 20, a, "three-holes.eps");

	return 0;
}

