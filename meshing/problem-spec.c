#include <stdio.h>
#include <math.h>
#include "array.h"
#include "problem-spec.h"

/*
double square_f(double x, double y)
{
	return 36*x*y*(1-x)*(1-y);
}

//Strength 9 Appears to be necesary for an exact solution
double square_f(double x, double y)
{
	return (40/3)*x*(1-pow(x, 2))*y*(1-pow(y, 3));
} */

double square_f(double x, double y)
{
	double Pi = 4.0*atan(1.0);
	return ((Pi*Pi)/4)*sin(Pi*x)*sin(Pi*y);
}

struct problem_spec *square(void)
{
	static struct problem_spec_point points[] = {
		{0,	0.0,	0.0,	FEM_BC_DIRICHLET},
		{1,	1.0,	0.0,	FEM_BC_DIRICHLET},
		{2,	1.0,	1.0,	FEM_BC_DIRICHLET},
		{3, 0.0,	1.0,	FEM_BC_DIRICHLET},
	};

	static struct problem_spec_segment segments[] = {
		{0,	0,	1,	FEM_BC_DIRICHLET},
		{1,	1,	2,	FEM_BC_DIRICHLET},
		{2,	2,	3,	FEM_BC_DIRICHLET},
		{3,	3,	0,	FEM_BC_DIRICHLET},
	};

	static struct problem_spec spec = {
		.points 	= points,
		.segments	= segments,
		.holes		= NULL,
		.npoints	= (sizeof points) / (sizeof points[0]),
		.nsegments 	= (sizeof segments) / (sizeof segments[0]),
		.nholes		= 0,
		.f		= square_f,
		.g		= NULL,
		.h		= NULL,
		.eta		= NULL,
		.u_exact	= NULL,
	};

	printf("Square\n");
	return &spec;
}

double three_holes_f(double x, double y)
{
	return pow(x, 2) + pow(y, 2);
}

struct problem_spec *three_holes(int n)
{
	double Pi = 4.0*atan(1.0);
	double r = 0.1;
	double len = 1;
	int k = 6;

	struct problem_spec_segment *segments;
	struct problem_spec_point *points;
	make_vector(segments, 3*n+k);
	make_vector(points, 3*n+k);

	struct problem_spec_hole *holes;
	make_vector(holes, 3);
	holes[0].x = -0.25;
	holes[0].y = -0.25;
	holes[1].x = -0.75;
	holes[1].y = -0.25;
	holes[2].x = -0.75;
	holes[2].y = -0.75;

	//Make the points and segments for the shape outline
	points[0] = (struct problem_spec_point) {0, 0.0,	0.0,	FEM_BC_DIRICHLET};
	points[1] = (struct problem_spec_point) {1, -len,	0.0,	FEM_BC_DIRICHLET};
	points[2] = (struct problem_spec_point) {2, -len,	-len,	FEM_BC_DIRICHLET};
	points[3] = (struct problem_spec_point) {3, -len/2,	-len,	FEM_BC_DIRICHLET};
	points[4] = (struct problem_spec_point) {4, -len/2,	-len/2,	FEM_BC_DIRICHLET};
	points[5] = (struct problem_spec_point) {5, 0.0,	-len/2,	FEM_BC_DIRICHLET};

	segments[0] = (struct problem_spec_segment) {0,	0,	1,	FEM_BC_DIRICHLET};
	segments[1] = (struct problem_spec_segment) {1,	1,	2,	FEM_BC_DIRICHLET};
	segments[2] = (struct problem_spec_segment) {2,	2,	3,	FEM_BC_DIRICHLET};
	segments[3] = (struct problem_spec_segment) {3,	3,	4,	FEM_BC_DIRICHLET};
	segments[4] = (struct problem_spec_segment) {4,	4,	5,	FEM_BC_DIRICHLET};
	segments[5] = (struct problem_spec_segment) {5,	5,	0,	FEM_BC_DIRICHLET};
	
	//Make the first hole at (-0.25, -0.25)
	for (int i = 0; i < n; i++, k++) {
		points[k].segment_no = k;
		points[k].x = r*cos(2*Pi/n*i) - holes[0].x;
		points[k].y = r*sin(2*Pi/n*i) - holes[0].y;
		points[k].bc = FEM_BC_DIRICHLET;

		segments[k].segment_no = k;
		segments[k].point_no_1 = k;
		segments[k].point_no_2 = k+1;
		segments[k].bc = FEM_BC_DIRICHLET;
	}
	segments[k-1].point_no_2 = k-n;

	//Make the second hole at (-0.75, -0.25)
	for (int i = 0; i < n; i++, k++) {
		points[k].segment_no = k;
		points[k].x = r*cos(2*Pi/n*i) - holes[1].x;
		points[k].y = r*sin(2*Pi/n*i) - holes[1].y;
		points[k].bc = FEM_BC_DIRICHLET;

		segments[k].segment_no = k;
		segments[k].point_no_1 = k;
		segments[k].point_no_2 = k+1;
		segments[k].bc = FEM_BC_DIRICHLET;
	}
	segments[k-1].point_no_2 = k-n;

	//Make the third hole at (-0.75, -0.75)
	for (int i = 0; i < n; i++, k++) {
		points[k].segment_no = k;
		points[k].x = r*cos(2*Pi/n*i) - holes[2].x;
		points[k].y = r*sin(2*Pi/n*i) - holes[2].y;
		points[k].bc = FEM_BC_DIRICHLET;

		segments[k].segment_no = k;
		segments[k].point_no_1 = k;
		segments[k].point_no_2 = k+1;
		segments[k].bc = FEM_BC_DIRICHLET;
	}
	segments[k-1].point_no_2 = k-n;

	struct problem_spec *spec;
	spec = xmalloc(sizeof *spec);

	spec -> points 		= points;
	spec -> segments	= segments;
	spec -> holes 		= holes;
	spec -> npoints 	= 3*n+6;
	spec -> nsegments	= 3*n+6;
	spec -> nholes		= 3;
	spec -> f 			= three_holes_f;
	spec -> g 			= NULL;
	spec -> eta 		= NULL;
	spec -> u_exact 	= NULL;

	printf("Square with holes\n");
	return spec;
}

/* ----------------------------Triangle With Hole-----------------------------*/

double f_triangle_with_hole(double x, double y)
{
	return 2.0/5*x*y*(2-x)*(2-y);
}

struct problem_spec *triangle_with_hole(void)
{
	static struct problem_spec_point points[] = {
		{0,	-1.0,	0.0,	FEM_BC_DIRICHLET},
		{1,	1.0,	0.0,	FEM_BC_DIRICHLET},
		{2,	0.0,	2.0,	FEM_BC_DIRICHLET},
		{3,	-0.25,	0.25,	FEM_BC_DIRICHLET},
		{4,	0.25,	0.25,	FEM_BC_DIRICHLET},
		{5,	0.25,	1.0,	FEM_BC_DIRICHLET},
		{6,	-0.25,	1.0,	FEM_BC_DIRICHLET},
	};

	static struct problem_spec_segment segments[] = {
		{0,	0,	1,	FEM_BC_DIRICHLET},
		{1,	1,	2,	FEM_BC_DIRICHLET},
		{2,	2,	0,	FEM_BC_DIRICHLET},
		{3,	3,	4,	FEM_BC_DIRICHLET},
		{4,	4,	5,	FEM_BC_DIRICHLET},
		{5,	5,	6,	FEM_BC_DIRICHLET},
		{6,	6,	3,	FEM_BC_DIRICHLET},
	};

	static struct problem_spec_hole holes[] = {
		{0.0,	0.75},
	};

	static struct problem_spec spec = {
		.points		= points,
		.segments	= segments,
		.holes		= holes,
		.npoints	= (sizeof points) / (sizeof points[0]),
		.nsegments	= (sizeof segments) / (sizeof segments[0]),
		.nholes		= (sizeof holes) / (sizeof holes[0]),
		.f		= *f_triangle_with_hole,
		.g		= NULL,
		.h		= NULL,
		.eta		= NULL,
		.u_exact	= NULL,
	};

	printf("Triangle with a hole\n");
	return &spec;
}

struct problem_spec *annulus(int n)
{
	struct problem_spec_point *points;
	double Pi = 4.0*atan(1.0);
	double a = 0.325, b = 2*a;

	make_vector(points, 2*n);

	//Make the inner circle points
	for (int i = 0; i < n; i++) {
		points[i].segment_no = i;
		points[i].x = a*cos(2*Pi/n*i);
		points[i].y = a*sin(2*Pi/n*i);
		points[i].bc = FEM_BC_DIRICHLET;
	}
	
	//Make the outer circle points
	for (int i = 0; i < n; i++) {
		points[n+i].segment_no = n + i;
		points[n+i].x = b*cos(2*Pi/n*i);
		points[n+i].y = b*sin(2*Pi/n*i);
		points[n+i].bc = FEM_BC_DIRICHLET;
	}

	//Make the circle segments
	struct problem_spec_segment *segments;
	make_vector(segments, 2*n);

	for (int i = 0; i < n; i++) {
		segments[i].segment_no = i;
		segments[i].point_no_1 = i;
		segments[i].point_no_2 = i+1;
		segments[i].bc = FEM_BC_DIRICHLET;
	}

	for (int i = 0; i < n; i++) {
		segments[n+i].segment_no = n+i;
		segments[n+i].point_no_1 = n+i;
		segments[n+i].point_no_2 = n+i+1;
		segments[n+i].bc = FEM_BC_DIRICHLET;
	}
	
	segments[n-1].point_no_2 -= n;
	segments[2*n-1].point_no_2 -= n;
	
	struct problem_spec_hole *holes;
	make_vector(holes, 1);
	holes[0].x = 0.0;
	holes[0].y = 0.0;

	struct problem_spec *spec;
	spec = xmalloc(sizeof *spec);

	spec -> points 		= points;
	spec -> segments	= segments;
	spec -> holes 		= holes;
	spec -> npoints 	= 2*n;
	spec -> nsegments	= 2*n;
	spec -> nholes		= 1;
	spec -> f 		= NULL;
	spec -> g 		= NULL;
	spec -> eta 		= NULL;
	spec -> u_exact 	= NULL;

	printf("An %d-sided annulus\n", n);
	return spec;
}

void free_spec(struct problem_spec *spec)
{
	if (spec != NULL) {
		free_vector(spec->points);
		free_vector(spec->segments);
		free_vector(spec->holes);
		free(spec);
	}
}
