RM      = /bin/rm -f
PGA_DIR = ${CURDIR}

# Defaults for various compiler settings.
# We leave CC and FC settings to the make defaults.
CPPFLAGS = -DFORTRANUNDERSCORE -D_REENTRANT
CFLAGS = -fPIC -Wall
FFLAGS =
LDFLAGS =
LIBS = -lm
INCLUDES =
RANLIB = ranlib
OPT = -O3
# Set this to no (or something different from yes) if you don't want
# shared libraries
SHAREDLIBS=yes

# A note on word-size: The old "WL" macro is the number of bits in an
# unsigned long (used for bit-arrays). You only need to specify this if
# you have an unusual architecture where
# sizeof(unsigned long) * 8
# is not the number of bits in the unsigned long (i.e. if you have
# bytes with != 8 bits).
# CPPFLAGS += -DWL=42

ifeq (,${MPI})
    # Debian sets up all symlinks for the default mpi installed
    ifeq (,$(shell pkg-config --cflags mpi-c || true))
        MPI = serial
    else
        MPI = default
        CC = mpicc
        FC = mpif77
        MPI_INC := $(shell pkg-config --cflags mpi-c)
        MPI_LIB := $(shell pkg-config --libs mpi-c)
    endif
endif

# We determine MPI compile-time options from MPI setting
# You may chose to not supply an MPI argument but instead set relevant
# compiler flags like MPI_INC and MPI_LIB by hand.
ifeq (${MPI},openmpi)
    # Openmpi is default on debian. The real path is
    # architecture-specific and there is no distribution-independent way
    # to really find out the architecture used. So you may end up
    # specifying MPI_INC and MPI_LIB by hand.
    MPI_INC := -I /usr/include/mpi
    MPI_LIB := -lmpi_mpifh -lmpi
else ifeq (${MPI},mpich)
    MPI_INCLUDE:=$(shell ls -1d /usr/include/mpich /usr/include/*/mpich 2> /dev/null)
    CC = mpicc.mpich
    FC = mpif77.mpich
    MPI_INC := -I $(MPI_INCLUDE)
    MPI_LIB := -lmpich
else ifeq (${MPI},mpich2)
    CC = mpicc.mpich
    FC = mpif77.mpich
    MPI_INC := -I /usr/include/mpich
    MPI_LIB := -lmpich
else ifeq (${MPI},lam)
    CC = mpicc.lam
    FC = mpif77.lam
    MPI_INC := -I /usr/include/lam
    MPI_LIB := -llam
else ifeq (${MPI},serial)
    MPI_INC := -I ${PGA_DIR}/fakempi
    MPI_LIB :=
    CPPFLAGS += -DFAKE_MPI
    MPI_STUB = $(PGA_LIB_DIR)/mpi_stub.o
endif

ifeq (${DEBUG},1)
    CFLAGS += -g
    PGA_LIB = pgapack-debug
else
    CFLAGS += ${OPT}
    PGA_LIB = pgapack
endif

# We still support the old "ARCH_TYPE" setting.
# For whatever it's worth these days. You probably want to leave these
# unspecified.
# Note: The old make framework set a architecture macro -D{ARCH_TYPE}.
# There is no architecture-specific code, so we leave this out.
ifeq (${ARCH_TYPE},sun4)
    CC = cc
    FC = f77
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},next)
    CC = cc
    FC =
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},rs6000)
    CC = cc
    FC = f77
else ifeq (${ARCH_TYPE},freebsd)
    CC = cc
    FC = f77
    FFLAGS = -w
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},irix)
    CC = cc
    FC = f77
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},alpha)
    CC = cc
    FC = f77
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},hpux)
    CC = cc
    FC = f77
else ifeq (${ARCH_TYPE},linux)
    CC = cc
    FC = f77
    FFLAGS = -w
    LDFLAGS += -s
else ifeq (${ARCH_TYPE},t3d)
    CC = /mpp/bin/cc
    FC = /mpp/bin/cf77
    CFLAGS = -T cray-t3d
    FFLAGS = -C cray-t3d -dp
    CPPFLAGS += -DFORTRANCAP
else ifeq (${ARCH_TYPE},powerchallenge)
    CC = cc
    FC = f77
    CFLAGS = -mips4 -fullwarn -64
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},challenge)
    CC = cc
    FC = f77
    CFLAGS = -mips2 -fullwarn -32
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},paragon)
    CC = cc
    FC = f77
    CFLAGS += -nx
    FFLAGS += -nx
    CPPFLAGS += -DFORTRANUNDERSCORE
else ifeq (${ARCH_TYPE},sp2)
    CC = cc
    FC = f77
else ifeq (${ARCH_TYPE},exemplar)
    CC = cc
    FC = fort77
    FFLAGS += -L/usr/lib -lU77
endif

ifeq (,${ARCH_TYPE})
    ARCH_TYPE = generic
endif

PGA_LIB_DIR = ${PGA_DIR}/lib/${ARCH_TYPE}-${MPI}
INCLUDES += -I ${PGA_DIR}/include ${MPI_INC}
LDFLAGS += -L ${PGA_LIB_DIR} -l ${PGA_LIB} ${MPI_LIB} ${LIBS}
CFLAGS += ${CPPFLAGS} ${INCLUDES}
FFLAGS += ${INCLUDES}
# Pass to sub-make
export CFLAGS
export LDFLAGS
export FFLAGS
export CC
export FC
export PGA_LIB
export MPI_STUB
export PGA_LIB_DIR
export RANLIB
export MPI
export SHAREDLIBS

all:
	mkdir -p $(PGA_LIB_DIR)
	$(MAKE) -C source
	$(MAKE) -C examples
	$(MAKE) -C test
	$(MAKE) -C docs

clean:
	$(MAKE) -C source   clean
	$(MAKE) -C examples clean
	$(MAKE) -C test     clean
	$(MAKE) -C docs     clean

clobber: clean
	rm -rf ${PGA_DIR}/lib

.PHONY: all clean clobber
