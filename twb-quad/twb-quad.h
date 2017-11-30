/* twb-quad.h
 *
 * Tables of points and weights for quadrature on a triangle.
 * 
 * From:
 * Taylor, Wingate and Bos, "Several new quadrature formulas
 * for polynomial integration in the triangle", Jan. 2005,
 * available at:
 *     http://arxiv.org/abs/math/0501496v2
 * A Feb. 2007 version updates a reference.  The tabular
 * data remains unchanged.
 *
 * Program by: Rouben Rostamian <rostamian@umbc.edu>
 * November 2005
 * Revised for Math 625: April 2006
 * Removed preprocessor macros: July 2007
 * Revised fall 2007: Added lambda3 in the qdat structure
 *
 * ------------------------------------------
 *
 * For integration on a triangle, a "quadrature table" is
 * a list of coordinates and weights such that the integral
 * of a function on the triangle may be approximated
 * by the weighted sum of the function values.
 *
 * Quadrature tables come in varying sizes.  Taylor, Wingate,
 * Bos (TWB) provide 14 quadrature tables. the smallest tables
 * consists of three points (and the corresponding weights).
 * The largest table consists of 120 points.  To each table
 * there corresponds a number d, called the quadrature's "strength".
 * The integration error is zero (within the floating point accuracy)
 * for all polynomials of degree up to d.  In TWB's tables, the strength
 * ranges from 2 to 25.
 *
 * To each quadrature table I have appended a dummy "point" of
 * weight -1 to serve as a sentinel to mark the end of the table.
 * The sentinel is not included in the count of integration points.
 *
 * The area of the TWB triangle is 2, defined as TWB_STANDARD_AREA below.
 * Scale as needed.
*/

#ifndef H_TWB_QUAD_H
#define H_TWB_QUAD_H

#define TWB_STANDARD_AREA 2.0

/* A quadrature table is an array of TWB_qdat structures. */
struct TWB_qdat {
	double lambda1;
	double lambda2;
	double weight;
	double lambda3;
};

/* twb_quad()

   Selects and returns a pointer to one of the 14 quadrature tables.

   Arguments:

        d [In/Out]

	   On entry, *d holds the desired quadrature strength.
	   If it is not equal of one of the 14 available
	   strengths, the table with the next largest strength
	   is selected.  If *d is greater than the largest
	   strength, then the table with the largest strength
	   is selected.  On return, *d is set to the strength
	   of the selected table.

       n [Out]

          The number of quadrature points in the selected table.

   Return value:

       Returns the selected table.

*/
struct TWB_qdat *twb_qdat(int *d, int *n);

#endif /* H_TWB_QUAD_H */

