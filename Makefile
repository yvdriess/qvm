CC = gcc

SOURCES = qvm.c

TARGETS = qvm

VPATH = sexp/lib
INCPATH = -I./sexp/include -I./
LIBPATH = -L./sexp/lib
LIBS = -lsexp -lquantum

CFLAGS = -Wall -O2 $(INCPATH) $(LIBPATH)

DEST_OBJS=$(SOURCES:.c=.o)

all:  qvm

qvm: $(DEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $(DEST_OBJS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGETS) $(DEST_OBJS)
