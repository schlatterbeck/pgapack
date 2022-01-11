TEST := $(shell echo *.data)
VRFY := $(patsubst %.data,%.vrfy,$(TEST))
PFX  := ../examples
PROG := $(PFX)/mgh/testprog $(PFX)/mgh/testprogde \
    $(PFX)/c/classic $(PFX)/c/dejong $(PFX)/c/example $(PFX)/c/maxbit   \
    $(PFX)/c/maxchar $(PFX)/c/maxint $(PFX)/c/name $(PFX)/c/namefull    \
    $(PFX)/c/udtstr                                                     \
    $(PFX)/fortran/classic $(PFX)/fortran/dejong $(PFX)/fortran/example \
    $(PFX)/fortran/maxbit $(PFX)/fortran/maxchar $(PFX)/fortran/maxint  \
    $(PFX)/fortran/name $(PFX)/fortran/namefull

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

test: unit $(VRFY)

unit:
	make -C unit

instverf.vrfy: instverf instverf.data
	./runtest instverf > instverf.vrfy
	test $$(grep 'was successful' instverf.vrfy | wc -l) -eq 5

%.vrfy: %.data $(PROG)
	./runtest $< | diff $< - > $@

%.data:
	./runtest $@ > $@

testclean:
	rm -f *.vrfy

clean: testclean
	rm -f instverf instverf.o

.PHONY: all clean unit
