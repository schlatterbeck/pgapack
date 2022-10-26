TEST_F := $(shell echo mgh*.data fortran*.data)
TEST_C := $(shell echo c_*.data deb_*.data nsga*.data seq*.data inst*.data)
PFX  := ../examples
PROG_F := \
    $(PFX)/mgh/testprog $(PFX)/mgh/testprogde $(PFX)/mgh/testprogdertr  \
    $(PFX)/fortran/classic $(PFX)/fortran/dejong $(PFX)/fortran/example \
    $(PFX)/fortran/maxbit $(PFX)/fortran/maxchar $(PFX)/fortran/maxint  \
    $(PFX)/fortran/name $(PFX)/fortran/namefull
PROG_C := $(PFX)/c/udtstr                                                 \
    $(PFX)/c/classic $(PFX)/c/dejong $(PFX)/c/example $(PFX)/c/maxbit   \
    $(PFX)/c/maxchar $(PFX)/c/maxint $(PFX)/c/name $(PFX)/c/namefull    \
# If we don't have a fortran compiler use only C examples
ifeq (,${FC})
    PROG := $(PROG_C)
    TEST := $(TEST_C)
else
    PROG := $(PROG_C) $(PROG_F)
    TEST := $(TEST_C) $(TEST_F)
endif
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

all: test

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
	make -C unit clobber

clean: testclean
	rm -f instverf instverf.o

.PHONY: all clean unit
