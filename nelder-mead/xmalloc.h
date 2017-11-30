#include <stdlib.h>

#ifndef H_MALLOC_H
#define H_MALLOC_H

void *malloc_or_exit(size_t nbytes, const char *name, int line);
#define xmalloc(nbytes) malloc_or_exit(nbytes, __FILE__, __LINE__)

#endif /* H_MALLOC_H */
