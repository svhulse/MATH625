#ifndef H_PROBLEM_SPEC_H
#define H_PROBLEM_SPEC_H
#define FEM_BC_DIRICHLET	2
#define FEM_BC_NEUMANN		3

struct problem_spec_point {
	int segment_no;
	double x;
	double y;
	int bc;
};

struct problem_spec_segment {
	int segment_no;
	int point_no_1;
	int point_no_2;
	int bc;
};

struct problem_spec_hole {
	double x;
	double y;
};

struct problem_spec {
	struct problem_spec_point *points;
	struct problem_spec_segment *segments;
	struct problem_spec_hole *holes;
	int npoints;
	int nsegments;
	int nholes;
	double (*f)(double x, double y);
	double (*g)(double x, double y);
	double (*h)(double x, double y);
	double (*eta)(double x, double y);
	double (*u_exact)(double x, double y);
};
#endif /*H_PROBLEM_SPEC*/
