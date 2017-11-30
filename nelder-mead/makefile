CC = gcc
CFLAGS = -Wall -pedantic -std=c99

OFILES = demo-2D.o nelder-mead.o xmalloc.o
TARGET = demo-2D

$(TARGET): $(OFILES)
		$(CC) $(OFILES) -o $@

clean:
		rm -f $(OFILES) $(TARGET)

# below this is the output of "gcc -MM":
demo-2D.o: demo-2D.c nelder-mead.h
nelder-mead.o: nelder-mead.c array.h xmalloc.h nelder-mead.h
xmalloc.o: xmalloc.c xmalloc.h
