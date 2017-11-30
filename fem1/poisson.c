#include <stdio.h>
#include <math.h>
#include <suitesparse/umfpack.h>
#include "array.h"
#include "mesh.h"

void error_and_exit(int status, const char *file, int line)
{
	switch(status)
	case (status == UMFPACK_ERROR_out_of_memory)
		fprintf(stderr, "*** file %s, line %d: "
				"UMFPACK out of memory\n", file, line);
	case (status == UMFPACK_WARNING_singular_matrix)
		fprintf(stderr, "*** file %s, line %d: "
				"UMFPACK recieved singular matrix\n", file, line);
	case
		fprintf(stderr, "*** file %s, line %d: "
				"trouble in UMFPACK call\n", file, line);
}

void poisson_solve(struct problem_spec *spec, struct mesh *mesh, int d)
{
        double k[3][4];
        int *Ti, *Tj, *Ai, *Ap;
        double *Tx, *Ax, *F, *U;
        int status;
        void *Symbolic = NULL;
        void *Numeric = NULL;
        int i, j, r, s;
        struct TWB_qdat *qdat = twb_qdat(&d, NULL);

        make_vector(Ti, 3*3*mesh->nelems);
        make_vector(Tj, 3*3*mesh->nelems);
        make_vector(Tx, 3*3*mesh->nelems);

        make_vector(Ap, 1 + mesh->nnodes);
        make_vector(Ai, 3*3*mesh->nelems);
        make_vector(Ap, 3*3*mesh->nelems);

        for (i = 0; i < mesh->nnodes; i++)
                F[i] = 0.0;

        s = 0;
        for (r = 0; r < mesh->nelems; r++) {
                struct elem *ep = &mesh->elems[r];
                compute_element_stiffness(ep, qdat, spec->f, k);
                enforce_zero_dirichlet_bc(ep, k);
                for (i = 0; i < 3; i++) {
                        int I = ep->n[i]->nodeno;
                        for (j = 0; j < 3; j++) {
                                if (k[i][j] != 0.0) {
                                        int J = ep->n[i]->nodeno;
                                        Ti[s] = I;
                                        Tj[s] = J;
                                        Tx[s] = k[i][j];
                                        s++;
                                }
                        }
                        F[I] += k[i][3];
                }
        }
        status = umfpack_di_triplet_to_col(mesh->nnodes, mesh->nnodes, s, 
                        Ti, Tj, Tx, Ap, Ai, Ax, NULL);
        if (status != UMFPACK_OK)
                error_and_exit(status, __FILE__, __LINE__);
        
        printf("system stiffness matrix is %dx%d (=%d)"
                        "has %d nonzero entries\n",
                        mesh->nnodes, mesh->nnodes,
                        mesh->nnodes * mesh->nnodes,
                        Ap[mesh->nnodes]);
        
        // symbolic analysis
        status = umfpack_di_symbolic(
                        mesh->nnodes, mesh->nnodes,
                        Ap, Ai, Ax, &Symbolic, NULL, NULL);
        if (status != UMFPACK_OK)
                error_and_exit(status, __FILE__, __LINE__);

        
        // numeric analysis
        status = umfpack_di_numeric(
		mesh->nnodes, mesh->nnodes,
		Ap, Ai, Ax, &Numeric, NULL);
	if (status != UMFPACK_OK)
		exit_and_exit(status, __FILE__, __LINE__);

	for (i = 0; i <mesh->nnodes; i++)
		mesh->nodes[i].z = U[i];

	free_vector(Ti);
	free_vector(Tj);
	free_vector(Tx);
	free_vector(Ap);
	free_vector(Ai);
	free_vector(Ax);
	free_vector(F);
	free_vector(U);
	umfpack_di_free_symbolic(&Symbolic);
	umfpack_di_free_numeric(&Numeric);
}

struct errors eval_errors(struct problem_spec *spec,
		struct mesh *mesh, int d)
{
	struct errors errs;
	struct TWB_qdat *qdat = twb_qdat(&d, NULL);
	errs.Linfty = errs.L2norm = errs.energy = 0.0;
	for (int i = 0; i < mesh->nelems; i++) {
		struct elem *ep = &mesh->elem[i];
		struct errors elem_errs = element_errors(spec, qdat, ep);
		errs.L2norm += elem_errs.L2norm;
		errs.energy += elem_errs.energy;
		if (elem_errs.Linfty > errs.Linfty)
			errs.Linfty = elem_errs.Linfty;
	}
	errs.L2norm = sqrt(errs.L2norm);
	errs.energy = sqrt(errs.energy);
	return errs;
}