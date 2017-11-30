CC = gcc
CFLAGS = -Wall -pedantic -std=c99
OFILES = mesh.o problem-spec.o twb-quad-demo.o twb-quad.o xmalloc.o \
	 plot-with-geomview.o
TARGET = twb-quad-demo

$(TARGET): $(OFILES) triangle.o
	$(CC) $(OFILES) triangle.o -o $@ -lm

clean:
	rm -f $(OFILES) $(TARGET)

#below this is the output of "gcc -MM":
twb-quad-demo.o: twb-quad-demo.c twb-quad.h problem-spec.h triangle.h
