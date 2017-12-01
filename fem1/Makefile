CC = gcc
CFLAGS = -Wall -pedantic -std=c99
OFILES = fem-demo.o poisson.o plot-with-geomview.o mesh.o problem-spec.o twb-quad.o xmalloc.o \
	 triangle.o
TARGET = fem-demo

$(TARGET): $(OFILES) triangle.o
	$(CC) $(OFILES) triangle.o -o $@ -lm -lumfpack

clean:
	rm -f $(OFILES) $(TARGET)

#below this is the output of "gcc -MM":
