OBJS = pop_test.o linalg_test.o integer_test.o utility_test.o
BIN  = $(patsubst %.o,%,$(OBJS))
DATA = $(patsubst %.o,%.data,$(OBJS))
VRFY = $(patsubst %.data,%.vrfy,$(DATA))
OUT  = $(patsubst %.data,%.out,$(DATA))
GCNO = $(DATA:.data=.gcno) $(DATA:_test.data=.gcno)
GCDA = $(DATA:.data=.gcda) $(DATA:_test.data=.gcda)
ROBJ = $(DATA:_test.data=.o)
ifeq (${CFLAGS},)
CFLAGS=-I../../include -I../../fakempi -g
endif
ifeq (${MPI},)
MPI=serial
endif

all: $(VRFY)

%_test: %_test.o %.o
	$(CC) -o $@ $^ ../../lib/generic-$(MPI)/*.a $(LDFLAGS) $(LIBS)

%.vrfy: % %.data
	$< > $*.out
	diff $*.out $*.data > $@

%.data:
	$* > $@

%.o: ../../source/%.c
	$(CC) $(CFLAGS) -DSTATIC= -c $< -o $@

clean: 
	rm -f $(OBJS) $(ROBJ) $(BIN) $(GCNO) $(GCDA)

clobber: clean
	rm -f $(VRFY) $(OUT)

# vim: set noet :
