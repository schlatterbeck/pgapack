OBJS = pop_test.o linalg_test.o
BIN  = $(patsubst %.o,%,$(OBJS))
DATA = $(patsubst %.o,%.data,$(OBJS))
VRFY = $(patsubst %.data,%.vrfy,$(DATA))
CFLAGS=-I../../include -I../../fakempi -g

all: $(VRFY)

%_test: %_test.o ../../source/%.c
	$(CC) $(CFLAGS) -DSTATIC= -o $@ $^ ../../lib/generic-serial/*.a -lm

%.vrfy: % %.data
	$< | diff $*.data - > $@

%.data:
	$* > $@

clean: 
	rm -f $(OBJS) $(BIN)

clobber: clean
	rm -f $(VRFY)
