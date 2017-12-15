#include <stdio.h>
#include <math.h>
#include "problem-spec.h"

static double heat1_exact(double x, double t)
{
	double Pi = 4.0*atan(1.0);
	return exp(-Pi*Pi/4*t) * cos(Pi/2*x);
}

static double heat1_ic(double x)
{
	return heat1_exact(x, 0);
}

static double heat1_bcL(double t)
{
	return heat1_exact(-1, t);
}

static double heat1_bcR(double t)
{
	return heat1_exact(1, t);
}

struct problem_spec *heat1(void)
{
	static struct problem_spec spec = {
		.a		= -1.0,
		.b		= 1.0,
		.ic 		= heat1_ic,
		.bcL 		= heat1_bcL,
		.bcR	 	= heat1_bcR,
		.u_exact 	= heat1_exact,
	};
	printf("problem heat1:\n");
	return &spec;
}

static double heat2_exact(double x, double t)
{
	double u;
	if (t > 0)
		u = erf(x/(2*sqrt(t)));
	else {
		if (x < 0.0)
			u = -1.0;
		else if (x > 0)
			u = +1.0;
		else
			u = 0.0;
	}
	return 0.4*u;
}

static double heat2_ic(double x)
{
	return heat2_exact(x, 0);
}

static double heat2_bcL(double t)
{
	return heat2_exact(-1, t);
}

static double heat2_bcR(double t)
{
	return heat2_exact(1, t);
}

struct problem_spec *heat2(void)
{
	static struct problem_spec spec = {
		.a		= -1.0,
		.b		= 1.0,
		.ic 		= heat2_ic,
		.bcL 		= heat2_bcL,
		.bcR	 	= heat2_bcR,
		.u_exact 	= heat2_exact,
	};
	printf("problem heat2:\n");
	return &spec;	
}

static double heat3_ic(double x)
{
	return fabs(x) < 0.4 ? 1.0 : 0.0;
}