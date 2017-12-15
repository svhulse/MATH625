#ifndef H_ARRAY_H
#define H_ARRAY_H

#include "xmalloc.h"

#define free_vector(x) do {free(x); x = NULL;} while (0)
#define make_vector(v, n)  ((v) = xmalloc((n) * sizeof *(v)))

#define make_matrix(a, m, n) do {           \
  make_vector(a, (m) + 1); 									\
  for (size_t make_matrix_i = 0; 						\
		make_matrix_i < (m);								    \
		make_matrix_i++)								        \
	make_vector((a)[make_matrix_i], (n));			\
  (a)[m] = NULL; 										        \
} while (0)

#define free_matrix(a) do { 								\
  if ((a) != NULL) { 										    \
    for (size_t free_matrix_i = 0;					\
		(a)[free_matrix_i] != NULL;	 						\
		free_matrix_i++)								        \
      free_vector((a)[free_matrix_i]); 			\
    free_vector(a); 										    \
    (a) = NULL; 										        \
  } 												                \
} while (0)

#define print_vector(fmt, v, n) do {								\
  for (int print_vector_i = 0; print_vector_i < (n); print_vector_i++) {			\
    printf(fmt, (v)[print_vector_i]); 								\
  } 												\
  printf("\n"); 										\
} while (0)

#define print_matrix(fmt, A, m, n) do {								\
  for (int print_matrix_i = 0; print_matrix_i < (m); print_matrix_i++) {			\
    print_vector(fmt, (A)[print_matrix_i], (n)); 						\
    printf("\n"); 										\
  } 												\
} while (0)

#endif /* H_ARRAY_H */
