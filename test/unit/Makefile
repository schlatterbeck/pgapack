OBJS = pop_test.o linalg_test.o integer_test.o utility_test.o
BIN  = $(patsubst %.o,%,$(OBJS))
DATA = $(patsubst %.o,%.data,$(OBJS))
VRFY = $(patsubst %.data,%.vrfy,$(DATA))
OUT  = $(patsubst %.data,%.out,$(DATA))
ifeq (${CFLAGS},)
CFLAGS=-I../../include -I../../fakempi -g
endif
ifeq (${MPI},)
MPI=serial
endif

all: $(VRFY)

%_test: %_test.o ../../source/%.c
	$(CC) $(CFLAGS) -DSTATIC= -o $@ $^ ../../lib/generic-$(MPI)/*.a -lm

%.vrfy: % %.data
	$< > $*.out
	diff $*.out $*.data > $@

%.data:
	$* > $@

clean: 
	rm -f $(OBJS) $(BIN)

clobber: clean
	rm -f $(VRFY) $(OUT)

# vim: set noet :
