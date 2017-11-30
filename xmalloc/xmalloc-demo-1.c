#include <stdio.h>
#include "xmalloc.h"

int main(void)
{
  double *x;
  x = malloc(1000);
  printf("Malloc of 1000 bytes executed\n");

  free(x);
  printf("%u bytes freed\n", (unsigned int)(sizeof(x)));
  x = NULL;

  x = malloc(0);
  printf("Malloc of 0 bytes executed\n");

  free(x);
  printf("Memmory freed\n");
  x = NULL;

  return EXIT_SUCCESS;
}
