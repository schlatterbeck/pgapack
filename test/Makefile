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
    $(PFX)/c/maxint-de
# If we don't have a fortran compiler use only C examples
ifeq (,${FC})
    PROG := $(PROG_C)
    TEST := $(TEST_C)
else
    PROG := $(PROG_C) $(PROG_F)
    TEST := $(TEST_C) $(TEST_F)
endif
VRFY := $(patsubst %.data,%.vrfy,$(TEST))
# Enable valgrind tests only for serial version
ifeq ($(MPI),serial)
MDBG := $(patsubst %.data,%.mdbg,$(TEST))
endif

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

test: unit $(VRFY) $(MDBG)

unit:
	make -C unit

instverf.vrfy: instverf instverf.data
	./runtest instverf > instverf.vrfy
	test $$(grep 'was successful' instverf.vrfy | wc -l) -eq 5

# Valgrind memory tests, we require no errors and no leaks
instverf.mdbg: instverf instverf.data
	./runtest -m instverf > instverf.vrfy-m 2> instverf.mdbg
	test $$(grep 'was successful' instverf.vrfy-m | wc -l) -eq 5
	grep -q 'ERROR SUMMARY: 0 errors from 0 contexts' instverf.mdbg
	grep -q 'All heap blocks were freed' instverf.mdbg

%.vrfy: %.data $(PROG)
	./runtest $< | diff $< - > $@

# Valgrind memory tests, we require no errors and no leaks
%.mdbg: %.data $(PROG)
	./runtest -m $< 2> $*.mdbg | diff $< - > $*.vrfy-m
	grep -q 'ERROR SUMMARY: 0 errors from 0 contexts' $*.mdbg
	grep -q 'All heap blocks were freed -- no leaks are possible' $*.mdbg

%.data:
	./runtest $@ > $@

testclean:
	rm -f *.vrfy *.vrfy-m *.mdbg
	make -C unit clobber

clean: testclean
	rm -f instverf instverf.o

.PHONY: all clean unit

# vim: set noet :
