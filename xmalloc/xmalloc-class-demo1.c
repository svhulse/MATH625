#include <stdio.h>
#include "xmalloc.h"

int main(void)
{
  double *x;
  int n = 100;

  printf("allocating a vector of &d double\n", n);

  x = xmalloc(n * sizeof *x);

  printf("freeing the allocated memory\n");
  free(x);

  return EXIT_SUCCESS;
}
