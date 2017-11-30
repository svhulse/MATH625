#include <stdio.h>
#include "array.h"
#include "nelder-mead.h"

#define REFLECT  1.0
#define EXPAND   2.0
#define CONTRACT 0.5
#define SHRINK   0.5

static inline void rank_vertices(double *y, int m, int *ia, int *iy, int *iz)
{
	double init_val_low = -999999;    //Must be lower than all values of the function

	double ia_val = y[0];
	double *iy_val = &init_val_low;
	double *iz_val = &init_val_low;

	*ia = 0;
	*iy = 0;
	*iz = 0;

	for (int i = 0; i < m; i++) {
		if (ia_val > y[i]) {
			ia_val = y[i];
			*ia = i;
		}

		if (*iz_val <= y[i]) {
			iy_val = iz_val;
			*iy = *iz;

			iz_val = &y[i];
			*iz = i;
		}

		else if (*iy_val <= y[i]) {
			iy_val = &y[i];
			*iy = i;
		}
	}
}

static void get_centroid(double **s, int n, int iz, double *C)
{
	//Reset the centroid values to 0
	for (int i = 0; i < n; i++)
		C[i] = 0;

	//Create a row vector of column sums excluding the iz row
	for (int i = 0; i < n+1; i++) {
		if (i != iz) {
			for (int j = 0; j < n; j++)
				C[j] = C[j] + s[i][j];
		}
	}

	//Divide each centroid value by n
	for (int i = 0; i < n; i++)
		C[i] = C[i] / n;
}

static inline void transform(double *P, double *Q, int n, double beta, double *R)
{
	for (int i = 0; i < n; i++)
		R[i] = (1 - beta)*P[i] + beta*Q[i];
}

static void shrink(double **s, int n, int ia)
{
	for (int i = 0; i < n+1; i++)
		if (i != ia)
			transform(s[ia], s[i], n, SHRINK, s[i]);
}

static inline void replace_row(double **s, int i, double **r)
{
	double *tmp = s[i];
	s[i] = *r;
	*r = tmp;
}

static int done(double **s, int n, double *y, int ia, int iz, double err2)
{
	double x_norm = 0;

	for (int i = 0; i < n; i++)
		x_norm = x_norm + (s[iz][i] - s[ia][i])*(s[iz][i] - s[ia][i]);

	if (x_norm < err2 && y[iz] - y[ia] < err2)
		return 1;
	else
		return 0;
}

int nelder_mead(struct nelder_mead *nm)
{
	double **s = nm->s;
	int n = nm->n;
	double h = nm->h;
	double tol = nm->tol;
	double err2 = (h*tol)*(h*tol);
	double *y, *C, *Pr, *Pe, *Pc;
	double yr, ye, yc;
	int ia, iy, iz;                     //indexies for best, worst, and second worst vertexies
	int simplex_to_be_freed = 0;
	int fevalcount;

	make_vector(y, n+1);     //vertex values
	make_vector(Pr, n);      //the reflected point x(r)
	make_vector(Pe, n);      //the expanded point x(e)
	make_vector(Pc, n);      //the contracted point x(c)
	make_vector(C, n);       //centoid of the face opposite vertex iz

	//Initialize vertex matrix is not already defined
	if (s == NULL) {
		make_matrix(s, n+1, n);
		simplex_to_be_freed = 1;

		for (int i = 0; i < n+1; i++)
			for (int j = 0; j < n; j++)
				s[i][j] = nm->x[j];

		for (int i = 1; i < n+1; i++)
			s[i][i-1] = s[i][i-1] + h;
	}

	//Calculate initial vertex values
	for(int i = 0; i < n+1; i++)
		y[i] = nm->f(s[i], n, nm->params);
	fevalcount = n+1;

	while (fevalcount <= nm->maxevals) {
		rank_vertices(y, n+1, &ia, &iy, &iz);

		//If covergence criteron met, store value and vector of minimized vertex
		if(done(s, n, y, ia, iz, err2)) {
			nm->minval = y[ia];
			
			for (int j = 0; j < n; j++)
				nm->x[j] = s[ia][j];
			break;
		}

		//Perform initial reflection, assign reflected iz value to yr
		get_centroid(s, n, iz, C);
		transform(C, s[iz], n, -REFLECT, Pr);
		yr = nm->f(Pr, n, nm->params);
		fevalcount++;

		if (yr < y[ia]) { //If the worst vertex becomes the best after reflection (Case 1)
			transform(C, Pr, n, EXPAND, Pe);
			ye = nm->f(Pe, n, nm->params);
			fevalcount++;

			if (ye < yr) {
				replace_row(s, iz, &Pe);
				y[iz] = ye;
			} else {
				replace_row(s, iz, &Pr);
				y[iz] = yr;
			}
		}

		else if (yr < y[iy]) { //If the worst vertex is better than the second best after reflection (Case 2)
			replace_row(s, iz, &Pr);
			y[iz] = yr;
		}

		else if (yr < y[iz]) {
			transform(C, Pr, n, CONTRACT, Pc);
			yc = nm->f(Pc, n, nm->params);
			fevalcount++;

			if (yc < yr) {
				replace_row(s, iz, &Pc);
				y[iz] = yc;
		       	} else {
				shrink(s, n, ia);
				for (int i = 0; i < n; i++)
					if (i != ia) {
						y[i] = nm->f(s[i], n, nm->params);
						fevalcount++;
					}
		       	}
		}

		else {
			transform(C, s[iz], n, CONTRACT, Pc);
			yc = nm->f(Pc, n, nm->params);
			fevalcount++;

			if (yc < y[iz]) {
				replace_row(s, iz, &Pc);
				y[iz] = yc;
		       	} else {
				shrink(s, n, ia);
				for (int i = 0; i < n; i++)
					if (i != ia) {
						y[i] = nm->f(s[i], n, nm->params);
						fevalcount++;
					}
			}
		}
	}

	//Free memory from each vector
	free_vector(y);
	free_vector(C);
	free_vector(Pr);
	free_vector(Pe);
	free_vector(Pc);

	if (simplex_to_be_freed)
		free_matrix(s);
	return fevalcount;
}
