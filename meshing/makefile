CC = gcc
CFLAGS = -Wall -pedantic -std=c99

OFILES = mesh-demo.o mesh-to-eps.o mesh.o problem-spec.o triangle.o xmalloc.o
TARGET = mesh-demo

$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@ -lm

clean:
	rm -f $(OFILES) $(TARGET)

#below this is the output of "gcc -MM":
mesh-demo.o: mesh-demo.c problem-spec.h mesh-to-eps.h mesh.h
