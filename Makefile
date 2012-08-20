CC = gcc

SOURCES = qvm.c

TARGETS = qvm

VPATH = sexp/lib
INCPATH = -I./sexp/include -I./
LIBPATH = #-L./sexp/lib
LIBS = -lsexp -lquantum
OFLAGS = -O2 -Wall #-O2
DFLAGS = # -g3
CFLAGS = $(OFLAGS) $(DFLAGS) $(INCPATH) $(LIBPATH) -std=c99

DEST_OBJS=$(SOURCES:.c=.o)

all:  qvm

qvm: $(DEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $(DEST_OBJS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGETS) $(DEST_OBJS)
