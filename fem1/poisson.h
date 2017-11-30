#ifndef H_POISSON_H
#define H_POISSON_H
#include "problem-spec.h"
#include "mesh.h"

struct errors {
        double Linfty;
        double L2norm;
        double energy;
};

void poisson_solve(struct problem_spec *spec,
                struct mesh *mesh, int d);

struct errors eval_errors(struct problem_spec *spec,
                struct mesh *mesh, int d);

#endif /* H_POISSON_H */