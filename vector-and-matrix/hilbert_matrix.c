#include <stdio.h>
#include "xmalloc.h"
#include "array.h"

double **hilbert_matrix(size_t n)
{
    double **H;

    make_matrix(H, n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        H[i][j] = 1.0 / (i + j + 1);

    return H;
}

int main(void)
{
  int n = 8;
  double **H = hilbert_matrix(n);

  print_matrix(H, n, n);
  free_matrix(H);

  return EXIT_SUCCESS;
}
