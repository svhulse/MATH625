OFILES = xmalloc.o plot-with-geomview.o
TARGET = plot-with-geomview
CFLAGS = -Wall -pedantic -std=c99 -O2
LIBS = -lm

$(TARGET): $(OFILES)
	cc $(OFILES) -o $@ $(LIBS)

clean:
	rm -f $(OFILES) $(TARGET)

plot-with-geomview.o: plot-with-geomview.c array.h xmalloc.h
