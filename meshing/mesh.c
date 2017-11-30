#include <stdio.h>
#include "triangle.h"
#include "xmalloc.h"
#include "array.h"
#include "mesh.h"
#include "problem-spec.h"

static struct triangulateio *problem_spec_to_triangle(struct problem_spec *spec)
{
	int i;
	struct triangulateio *in = xmalloc(sizeof *in);

	/* process points */
	in->numberofpoints = spec->npoints;
	make_vector(in->pointlist, 2 * in->numberofpoints);
	for (i = 0; i < in->numberofpoints; i++) {
		in->pointlist[2*i]   = spec->points[i].x;
		in->pointlist[2*i+1] = spec->points[i].y;
	}

	make_vector(in->pointmarkerlist, in->numberofpoints);
	for (i = 0; i < in->numberofpoints; i++)
		in->pointmarkerlist[i] = spec->points[i].bc;

	in->numberofpointattributes = 0;

	/* process segments */
	in->numberofsegments = spec->nsegments;
	make_vector(in->segmentlist, 2 * in->numberofsegments);
	for (i = 0; i < in->numberofsegments; i++) {
		in->segmentlist[2*i]   = spec->segments[i].point_no_1;
		in->segmentlist[2*i+1] = spec->segments[i].point_no_2;
	}

	make_vector(in->segmentmarkerlist, in->numberofsegments);
	for (i = 0; i < in->numberofsegments; i++)
		in->segmentmarkerlist[i] = spec->segments[i].bc;

	/* process holes */
	in->numberofholes = spec->nholes;
	if (in->numberofholes != 0) {
		make_vector(in->holelist, 2 * in->numberofholes);
		for (i = 0; i < in->numberofholes; i++) {
			in->holelist[2*i]   = spec->holes[i].x;
			in->holelist[2*i+1] = spec->holes[i].y;
		}
	} else
		in->holelist = NULL;

	in->numberofregions = 0;

	return in;
}

static struct triangulateio *do_triangulate(struct triangulateio *in, double a)
{
	char *opts_string = "Qzpeq30a%f";
	char opts[64];
	struct triangulateio *out = xmalloc(sizeof *out);

	out->pointlist = NULL;
	out->pointmarkerlist = NULL;
	out->edgelist = NULL;
        out->edgemarkerlist = NULL;
	out->trianglelist = NULL;
	out->segmentlist = NULL;
        out->segmentmarkerlist = NULL;
	sprintf(opts, opts_string, a);
	triangulate(opts, in, out, NULL);

	return out;
}

/* This is ugly; there ought to be a better way but I don't see it */
static void assign_elem_edges(
		struct elem *elems, int nelems,
		struct edge *edges, int nedges)
{
	for (int r = 0; r < nelems; r++) {
		for (int i = 0; i < 3; i++) {	/* i: vertex index */
			int j = (i+1)%3;
			int k = (i+2)%3;
			int n1 = elems[r].n[j]->nodeno;
			int n2 = elems[r].n[k]->nodeno;
			for (int s = 0; s < nedges; s++) {
				int m1 = edges[s].n[0]->nodeno;
				int m2 = edges[s].n[1]->nodeno;
				if ((m1 == n1 && m2 == n2)
						|| (m1 == n2 && m2 == n1)) {
					elems[r].e[i] = &edges[s];
					break;
				}
			}
		}
	}
}

static void set_element_edge_vectors(struct elem *ep)
{
	for (int i = 0; i < 3; i++) {	/* i: vertex index */
		int j = (i+1)%3;
		int k = (i+2)%3;
		ep->ex[i] = ep->n[k]->x - ep->n[j]->x;
		ep->ey[i] = ep->n[k]->y - ep->n[j]->y;
	}
}

static void set_element_area(struct elem *ep)
{
	ep->area = (ep->ex[0]*ep->ey[1] - ep->ex[1]*ep->ey[0])/2.0;
}

static void set_edge_vectors_and_areas(struct elem *elems, int nelems)
{
	for (int i = 0; i < nelems; i++) {
		struct elem *ep = &elems[i];
		set_element_edge_vectors(ep); 
		set_element_area(ep); 
	}
}

static struct mesh *triangle_to_mesh(struct triangulateio *out)
{
	struct node *nodes;
	struct edge *edges;
	struct elem *elems;
	int i, nnodes, nedges, nelems;
	struct mesh *mesh = xmalloc(sizeof *mesh);

	nnodes = out->numberofpoints;
	make_vector(nodes, nnodes);
	for (i = 0; i < out->numberofpoints; i++) {
		nodes[i].nodeno = i;
		nodes[i].x = out->pointlist[2*i];
		nodes[i].y = out->pointlist[2*i+1];
		nodes[i].z = 0.0;
		nodes[i].bc = out->pointmarkerlist[i];
	}

	nedges = out->numberofedges;
	make_vector(edges, nedges);
	for (i = 0; i < nedges; i++) {
		edges[i].edgeno = i;
		edges[i].n[0] = &nodes[out->edgelist[2*i]];
		edges[i].n[1] = &nodes[out->edgelist[2*i+1]];
		edges[i].bc = out->edgemarkerlist[i];
	}

	nelems = out->numberoftriangles;
	make_vector(elems, nelems);
	for (i = 0; i < nelems; i++) {
		elems[i].elemno = i;
		elems[i].n[0] = &nodes[out->trianglelist[3*i]];
		elems[i].n[1] = &nodes[out->trianglelist[3*i+1]];
		elems[i].n[2] = &nodes[out->trianglelist[3*i+2]];
	}

	assign_elem_edges(elems, nelems, edges, nedges);
	set_edge_vectors_and_areas(elems, nelems);

	mesh->nnodes = nnodes;
	mesh->nedges = nedges;
	mesh->nelems = nelems;
	mesh->nodes = nodes;
	mesh->edges = edges;
	mesh->elems = elems;

	return mesh;
}

static void free_triangle_in_structure(struct triangulateio *in)
{
	free_vector(in->pointlist);
	free_vector(in->pointmarkerlist);
	free_vector(in->segmentlist);
	free_vector(in->segmentmarkerlist);
	free_vector(in->holelist);
	free(in);
}

static void free_triangle_out_structure(struct triangulateio *out)
{

	free(out->pointlist);
	free(out->pointmarkerlist);
	free(out->edgelist);
	free(out->edgemarkerlist);
	free(out->trianglelist);
	free(out->segmentlist);
	free(out->segmentmarkerlist);
	free(out);
}

struct mesh *make_mesh(struct problem_spec *spec, double a)
{
	struct triangulateio *in, *out;
	struct mesh *mesh;

	in = problem_spec_to_triangle(spec);
	out = do_triangulate(in, a);
	mesh = triangle_to_mesh(out);
	free_triangle_in_structure(in);
	free_triangle_out_structure(out);
	return mesh;
}

void free_mesh(struct mesh *mesh)
{
	if (mesh == NULL)
		return;

	free_vector(mesh->nodes);
	free_vector(mesh->edges);
	free_vector(mesh->elems);
	free(mesh);
}

