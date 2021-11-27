TEST := $(shell echo *.data)
VRFY := $(patsubst %.data,%.vrfy,$(TEST))

all: instverf 

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

test: $(VRFY)

instverf.vrfy: instverf instverf.data
	./instverf > instverf.vrfy
	test $$(grep 'was successful' instverf.vrfy | wc -l) -eq 5

%.vrfy: %.data
	./runtest $< | diff - $< > $@

testclean:
	rm -f *.vrfy

clean: testclean
	rm -f instverf instverf.o

.PHONY: all clean
