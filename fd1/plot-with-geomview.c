#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"

double f(double x, double t)
{
	return (1-x*x)*exp(-t);
}

void plot_curve(FILE *fp, double *v, int n, int steps, int k)
{
	double t = (double)k / steps;
	for (int j = 0; j < n+2; j++) {
		double x = (double)j / (n+1);
		fprintf(fp, "%g %g %g\n", x, t, v[j]);
	}
}

int main(void)
{
	double a = -1.0, b = 1.0;
	double T = 1.0;
	int n = 15;
	int steps = 10;
	double *v;

	double dt = T / steps;
	double dx = (b-a) / (n+1);

	FILE *fp;
	fp = fopen("demo.gv", "w");
	if (fp == NULL) {
		fprintf(stderr, "Unable to open file for writing\n");
		exit(EXIT_FAILURE);
	}

	make_vector(v, n+2);

	fprintf(fp, "# my demo geomview file\n");
	fprintf(fp, "{ appearance { +edge } \n");
	fprintf(fp, "MESH %d %d\n", n+2, steps+1);

	for (int i = 0; i <= steps; i++) {
		double t = i * dt;
		for (int j = 0; j < n+2; j++) {
			double x = a + dx*j;
			v[j] = f(x, t);
		}
		plot_curve(fp, v, n, steps, i);
	}

	fprintf(fp, "}\n");
	fclose(fp);

	return EXIT_SUCCESS;
}