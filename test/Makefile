TEST := $(shell echo *.data)
VRFY := $(patsubst %.data,%.vrfy,$(TEST))

# You can override this on the command-line.
# Mainly needed for running tests in parallel mode with mpirun.
ifeq (,${MPI_MACHINEFILE})
    MPI_MACHINEFILE := $(HOME)/.mpi-$(MPI)-machinefile
endif
export MPI_MACHINEFILE
ifeq (,${MPI_NUM_PROCESS})
    MPI_NUM_PROCESS := 4
endif
export MPI_NUM_PROCESS

all: instverf 

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

test: $(VRFY)

instverf.vrfy: instverf instverf.data
	./runtest instverf > instverf.vrfy
	test $$(grep 'was successful' instverf.vrfy | wc -l) -eq 5

%.vrfy: %.data
	./runtest $< | diff $< - > $@

testclean:
	rm -f *.vrfy

clean: testclean
	rm -f instverf instverf.o

.PHONY: all clean
