/*
The "physical" coordinates range over xmin..xmax, and ymin..ymax.
The postscript bounding box is (0,0) to (W,H), where max(W,H) = D is
define below.

The physical domain is scaled to the postscript domain by a factor of s.
Thus, the width of the image area is s*w.  We add margins of p*s*w
from the left and right, so that 2*p*s*w + s*w = W.  Similarly,
2*p*s*h + s*h = H.  Therefore: W = (1+2*p)*s*w, and H = (1+2*p)*s*h, and
max(W,H) = (1+2*p)*s*max(w,h).  Consequently, D = (1+2*p)*s*d, where
d=max(w,h).  From here we get:
   s = D/((1+2*p)*d).
In particular, (1+2*p)s = D/d, therefore W = D/d*w and H = D/d*h.

The linear mapping from the physical domain (x,y) to the postscript
domain (X,Y) takes xmin to the images left margin, that is p*s*w,
therefore X = p*s*w + s*(x - xmin).  Similarly Y = p*s*h + s*(y - ymin).

These are all incorporated into the function mesh_to_eps().

2012-12-29
*/

#include <stdio.h>
#include <time.h>
#include "mesh-to-eps.h"

#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define D 400	/* the larger of the W & H of the bounding box */

void mesh_to_eps(struct mesh *mesh, char *outfile)
{
	FILE *fp;
	struct node *nodes = mesh->nodes;
	struct edge *edges = mesh->edges;
	struct elem *elems = mesh->elems;
	time_t now;
	double xmin, xmax, ymin, ymax, w, h, W, H, d, s;
	double p = 0.01;
	int i, k;

	if ((fp = fopen(outfile, "w")) == NULL) {
		fprintf(stderr, "cannot open file %s for writing\n", outfile);
		return;
	}

	xmin = xmax = nodes[0].x;
	ymin = ymax = nodes[0].y;
	for (i = 1; i < mesh->nnodes; i++) {
		if (nodes[i].x < xmin)
			xmin = nodes[i].x;
		else if(nodes[i].x > xmax)
			xmax = nodes[i].x;
		if (nodes[i].y < ymin)
			ymin = nodes[i].y;
		else if(nodes[i].y > ymax)
			ymax = nodes[i].y;
	}
	w = xmax - xmin;
	h = ymax - ymin;
	d = MAX(w, h);
	W = D/d*w;		/* bounding box width */
	H = D/d*h;		/* bounding box height */
	s = D/((1+2*p)*d);	/* scale */

	fputs("%!PS-Adobe-3.0 EPSF-3.0\n", fp);
	fprintf(fp, "%%%%BoundingBox: 0 0 %g %g\n", W, H);
	fprintf(fp, "%%%%Title: (%s)\n", outfile);
	fprintf(fp, "%%%%Creator: C Projects, %s\n", __FILE__);
	now = time(NULL);
	fprintf(fp, "%%%%CreationDate: %s", ctime(&now));
	fputs("%%EndComments\n", fp);

	fputs("gsave\n", fp);
	fputs("1 1 0 setrgbcolor\n", fp);
	for (i = 0; i < mesh->nelems; i++) {
		double x[3], y[3];
		for (k = 0; k < 3; k++) {
			x[k] = p*s*w + s*(elems[i].n[k]->x - xmin);
			y[k] = p*s*h + s*(elems[i].n[k]->y - ymin);
		}
		fprintf(fp, "%g %g moveto %g %g lineto %g %g lineto "
				"closepath fill\n",
				x[0], y[0],
				x[1], y[1],
				x[2], y[2]);
	}

	fputs("0 setgray\n", fp);
	for (i = 0; i < mesh->nedges; i++) {
		double x1 = edges[i].n[0]->x;
		double y1 = edges[i].n[0]->y;
		double x2 = edges[i].n[1]->x;
		double y2 = edges[i].n[1]->y;
		double X1 = p*s*w + s*(x1 - xmin);
		double Y1 = p*s*h + s*(y1 - ymin);
		double X2 = p*s*w + s*(x2 - xmin);
		double Y2 = p*s*h + s*(y2 - ymin);
		fprintf(fp, "%g %g moveto %g %g lineto stroke\n",
			X1, Y1, X2, Y2);
	}
	fputs("grestore\n", fp);
	fputs("showpage\n", fp);
	fputs("%%EOF\n", fp);
	fclose(fp);
	fprintf(stderr, "postscript output written to file %s\n", outfile);
}

