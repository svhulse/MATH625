#ifndef H_NELDER_MEAD_H
#define H_NELDER_MEAD_H

struct nelder_mead {
  double (*f)(double *x, int n, void *params);
  int n;
  double **s;
  double *x;
  double h;
  double tol;
  int maxevals;
  double minval;
  void *params;
};

int nelder_mead(struct nelder_mead *nm);
#endif /*H_NELDER_MEAD*/
