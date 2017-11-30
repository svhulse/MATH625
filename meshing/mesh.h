#ifndef H_MESH_H
#define H_MESH_H
#include "problem-spec.h"

struct node {
	int nodeno;
	double x;
	double y;
	double z;
	int bc;
};

struct edge {
	int edgeno;
	struct node *n[2];	//pointers to the edge's endpoints
	int bc;
};

struct elem {
	int elemno;
	struct node *n[3];	//pointers to the element's nodes
	struct edge *e[3];	//pointers to the element's edges
	double ex[3], ey[3];	//x and y components of the edge vectors
	double area;		//element's area
};

struct mesh {
	struct node *nodes;	//the array of node structures
	struct edge *edges;	//the array of edge structures
	struct elem *elems;	//the array of elem structures
	int nnodes;
	int nedges;
	int nelems;
};

struct mesh *make_mesh(struct problem_spec *spec, double a);
void free_mesh(struct mesh *mesh);
#endif /*H_MESH_H*/
