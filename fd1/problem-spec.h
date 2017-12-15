#ifndef H_PROBLEM_SPEC_H
#define H_PROBLEM_SPEC_H

struct problem_spec {
	double a;
	double b;
	double (*ic)(double x);
	double (*bcL)(double t);
	double (*bcR)(double t);
	double (*u_exact)(double x, double t);
};
#endif