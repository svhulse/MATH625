#ifndef H_ARRAY_H
#define H_ARRAY_H

#include "xmalloc.h"

#define free_vector(x) do {free(x); x = NULL;} while (0)
#define make_vector(v,n)  ((v) = xmalloc((n) * sizeof *(v)))

#define make_matrix(a, m, n) do { \
  make_vector(a, m + 1); \
  \
  for (size_t make_matrix_counter = 0; make_matrix_counter < m; make_matrix_counter++) \
    make_vector(a[make_matrix_counter], n); \
  \
  a[m] = NULL; \
} while (0)

#define free_matrix(a) do { \
  if (a != NULL) { \
    for (size_t free_matrix_counter = 0; a[free_matrix_counter] != NULL; free_matrix_counter++) \
      free_vector(a[free_matrix_counter]); \
    free_vector(a); \
    a = NULL; \
  } \
} while (0)

#define print_vector(A, n) do { \
  for (int print_vector_i = 0; print_vector_i < n; print_vector_i++) { \
    printf("%8.3f", A[print_vector_i]); \
  } \
  printf("\n"); \
} while (0)

#define print_matrix(A, m, n) do { \
  for (int print_matrix_i = 0; print_matrix_i < m; print_matrix_i++) { \
    for (int print_matrix_j = 0; print_matrix_j < n; print_matrix_j++) \
      printf("%8.3f", A[print_matrix_i][print_matrix_j]); \
    printf("\n"); \
  } \
} while (0)

#endif /* H_ARRAY_H */
