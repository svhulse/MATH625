/* plot-with-geomview.c
 *
 * Rouben Rostamian <rostamian@umbc.edu>
 * 2007-12-07
 *
 * Notes:  Within the geomview window press Pa to bring up the appearance
 * panel.  Experiment with various settings, e.g., turn the edges off then
 * select smooth shading.
*/
#include <stdio.h>
#include "plot-with-geomview.h"

#define REDHUE(s)       (huefunc((s)-23.0/64))
#define GREENHUE(s)     (huefunc((s) -7.0/64))
#define BLUEHUE(s)      (huefunc((s) +9.0/64))

/* This defines a continuous, piecewise-linear function defined on
   the entire real line.  Its graph is trapezoid-shaped and connects
   the points (0,0) (16/64,1), (33/64,1), (49/64,0).  The function
   is zero outside the support of the trapezoid.
   In combination with REDUE, GREENHUE, BLUEHUE macros defined above,
   produces an RGB scale used for our ZHUE rendering which coincides
   with Matlab's JET colormap.
*/
static float huefunc(float s)
{
	if (s < 0)
		return 0.0;
	else if (s < 16.0/64)
		return 4.0*s;
	else if (s < 33.0/64)
		return 1.0;
	else if (s < 49.0/64)
		return 1.0 - 4.0*(s-33.0/64);
	else
		return 0.0;
}

/* computes the min and max of the z values of the nodes of a mesh */
static void get_range(struct mesh *mesh, double *min, double *max)
{
	int i;

	*min = *max = mesh->nodes[0].z;
	for (i = 1; i < mesh->nnodes; i++) {
		double z = mesh->nodes[i].z;
		if (z < *min)
			*min = z;
		else if (z > *max)
			*max = z;
	}
}

/* In a COFF file, colors are specified at polygon vertices.
   These are interpolated smoothly to shade the polygon.
*/
static void draw_triangles_zhue(struct mesh *mesh, FILE *fp)
{
	double min, max, s;
	int i;

	get_range(mesh, &min, &max);

	fputs("COFF\n", fp);
	fprintf(fp, "%d %d %d\n", mesh->nnodes, mesh->nelems, mesh->nedges);

	for (i = 0; i < mesh->nnodes; i++) {
		struct node *np = &mesh->nodes[i];
		if (max - min < 1.0e-7)
			s = 0.5;
		else
			s = (np->z - min)/(max - min);	/* color scale */
		fprintf(fp, "%g %g %g %g %g %g 1.0\n",
				np->x, np->y, np->z,
				REDHUE(s), GREENHUE(s), BLUEHUE(s));
	}

	for (i = 0; i < mesh->nelems; i++) {
		struct elem *ep = &mesh->elems[i];
		fprintf(fp, "3 %d %d %d\n",
				ep->n[0]->nodeno,
				ep->n[1]->nodeno,
				ep->n[2]->nodeno);
	}
}

void plot_with_geomview_zhue(struct mesh *mesh, char *filename)
{
	FILE *fp;

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "unable to open file %s for writing\n",
				filename);
		return;
	}

	fputs("{ appearance { +edge +csmooth }\n", fp);
	draw_triangles_zhue(mesh, fp);
	fputs("}\n", fp);
	fclose(fp);
	printf("geomview output written to file %s\n", filename);
}

/* --- Monochrome plotting here --------------------------------- */

/* In an OFF file, the triangles may be specified as
   3 n1 n2 n3 R G B
   where n1 n2 n3 are the triangle's vertex numbers and R G B
   (each between 0 and 1) define the triangle's color.  Alternatively,
   we may omit the R G B, as in:
   3 n1 n2 n3
   in which case the triangle will be shown in a default silvery color.
*/
static void draw_triangles_mono(struct mesh *mesh, FILE *fp)
{
	int i;

	fputs("OFF\n", fp);
	fprintf(fp, "%d %d -1\n", mesh->nnodes, mesh->nelems);

	for (i = 0; i < mesh->nnodes; i++) {
		struct node *np = &mesh->nodes[i];
		fprintf(fp, "%g %g %g\n", np->x, np->y, np->z);
	}

	for (i = 0; i < mesh->nelems; i++) {
		struct elem *ep = &mesh->elems[i];
		fprintf(fp, "3 %d %d %d ",
				ep->n[0]->nodeno,
				ep->n[1]->nodeno,
				ep->n[2]->nodeno);
		fputs("1.0 1.0 0.0\n", fp);
	}
}

void plot_with_geomview_mono(struct mesh *mesh, char *filename)
{
	FILE *fp;

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "unable to open file %s for writing\n",
				filename);
		return;
	}

	fputs("{ appearance { +edge }\n", fp);
	draw_triangles_mono(mesh, fp);
	fputs("}\n", fp);
	fclose(fp);
	printf("geomview output written to file %s\n", filename);
}

