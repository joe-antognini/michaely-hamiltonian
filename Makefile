# installation prefix; executables are installed in $(PREFIX)/bin
PREFIX = $(HOME)

# test for ccache
CCACHE = $(shell which ccache 2>/dev/null)

ifneq ($(CCACHE),)
CC = ccache gcc
else
CC = gcc
endif

# test for architecture
UNAME = $(shell uname)

ifeq ($(UNAME),Linux)
CFLAGS = -Wall -O3
LIBFLAGS = -lgsl -lgslcblas -lm
else
ifeq ($(UNAME),Darwin)
CFLAGS = -Wall -O3 -I/sw/include -I/sw/include/gnugetopt -L/sw/lib
LIBFLAGS = -lgsl -lgslcblas -lgnugetopt -lm
else
CFLAGS = -Wall -O3
LIBFLAGS = -lgsl -lgslcblas -lm
endif
endif

# the core fewbody objects
FEWBODY_OBJS = fewbody.o fewbody_classify.o fewbody_coll.o fewbody_hier.o \
	fewbody_int.o fewbody_io.o fewbody_isolate.o fewbody_ks.o \
	fewbody_nonks.o fewbody_scat.o fewbody_utils.o 

all: binsingle 

binsingle: binsingle.o $(FEWBODY_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

binsingle.o: binsingle.c binsingle.h fewbody.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

install: binsingle
	mkdir -p $(PREFIX)/bin/
	install -m 0755 binsingle $(PREFIX)/bin/

uninstall:
	rm -f $(PREFIX)/bin/binsingle

clean:
	rm -f $(FEWBODY_OBJS) binsingle.o binsingle

mrproper: clean
	rm -f *~ *.bak *.dat ChangeLog
	rm -f */*~ */*.bak */*.dat
