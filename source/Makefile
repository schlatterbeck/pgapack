OBJS = $(PGA_LIB_DIR)/binary.o        \
      $(PGA_LIB_DIR)/char.o          \
      $(PGA_LIB_DIR)/cmdline.o       \
      $(PGA_LIB_DIR)/create.o        \
      $(PGA_LIB_DIR)/cross.o         \
      $(PGA_LIB_DIR)/debug.o         \
      $(PGA_LIB_DIR)/duplcate.o      \
      $(PGA_LIB_DIR)/evaluate.o      \
      $(PGA_LIB_DIR)/fitness.o       \
      $(PGA_LIB_DIR)/integer.o       \
      $(PGA_LIB_DIR)/linalg.o        \
      $(PGA_LIB_DIR)/mutation.o      \
      $(PGA_LIB_DIR)/parallel.o      \
      $(PGA_LIB_DIR)/pga.o           \
      $(PGA_LIB_DIR)/pop.o           \
      $(PGA_LIB_DIR)/random.o        \
      $(PGA_LIB_DIR)/real.o          \
      $(PGA_LIB_DIR)/report.o        \
      $(PGA_LIB_DIR)/restart.o       \
      $(PGA_LIB_DIR)/select.o        \
      $(PGA_LIB_DIR)/stop.o          \
      $(PGA_LIB_DIR)/system.o        \
      $(PGA_LIB_DIR)/user.o          \
      $(PGA_LIB_DIR)/utility.o $(PGA_LIB_DIR)/f2c.o $(MPI_STUB)

INSTALL=install
INSTALL_PLAIN=install -m 0644

########################################################
#                    The Target                        #
########################################################

v_git           := $(shell git describe)
MAJOR 		:= $(firstword $(subst ., ,$(v_git)))
MINOR 		:= $(word 2,$(subst ., ,$(subst -, ,$(v_git))))
PATCHLEVEL      := $(or $(word 2,$(subst -, ,$(v_git))),0)
VERSION         := $(MAJOR).$(MINOR).$(PATCHLEVEL)
ifeq (${PCFILE},)
PCFILE          := pgapack.pc
endif
# Target triplet
ifeq (${TARGET_TRIPLET},)
TARGET_TRIPLET:=$(shell ./triplet)
endif
REL_LIBDIR=lib/$(TARGET_TRIPLET)
INSTALL_LIBDIR=$(PREFIX)/$(REL_LIBDIR)
INSTALL_INCDIR=$(PREFIX)/include/pgapack-$(MPI)
ifeq (${SHAREDLIBS},yes)
SHLIBFILE = lib$(PGA_LIB).so.$(VERSION)
SHLIB = $(PGA_LIB_DIR)/$(SHLIBFILE)
SHARED = -shared -Wl,-soname,lib$(PGA_LIB).so.$(MAJOR)
SHLIB_INSTALL = $(INSTALL) $(SHLIB) $(INSTALL_LIBDIR) && \
    ln -sf $(SHLIBFILE) $(INSTALL_LIBDIR)/lib$(PGA_LIB).so
else
SHLIB_INSTALL = true
endif

ifeq (${MPI},serial)
FAKEMPI_OPT='-I$${includedir}/fakempi'
REQUIRES=
else
ifeq (${MPI},openmpi)
REQUIRES=Requires: ompi
else ifeq (${MPI},mpich)
REQUIRES=Requires: mpich
else ifeq (${MPI},lam)
# lam doesn't seem to have a .pc file
REQUIRES=
else
REQUIRES=Requires: mpi
endif
endif

all: $(PGA_LIB_DIR)/lib$(PGA_LIB).a $(SHLIB) $(PCFILE)

$(PGA_LIB_DIR)/lib$(PGA_LIB).a: $(OBJS)
	ar rv $(PGA_LIB_DIR)/lib$(PGA_LIB).a $(OBJS)
	$(RANLIB) $(PGA_LIB_DIR)/lib$(PGA_LIB).a

ifeq (${SHAREDLIBS},yes)
$(SHLIB): $(OBJS)
	$(CC) $(SHARED) $(OBJS) -o $(SHLIB)
	ln -s $(SHLIBFILE) $(PGA_LIB_DIR)/lib$(PGA_LIB).so.$(MAJOR)
	ln -s $(SHLIBFILE) $(PGA_LIB_DIR)/lib$(PGA_LIB).so
endif

$(PGA_LIB_DIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(PCFILE)

install: all
	$(INSTALL) -d $(INSTALL_LIBDIR)/pkgconfig
	$(INSTALL) -d $(INSTALL_INCDIR)
	$(SHLIB_INSTALL)
	$(INSTALL) -T $(PGA_LIB_DIR)/lib$(PGA_LIB).a \
	    $(INSTALL_LIBDIR)/lib$(PGA_LIB).a
	$(INSTALL_PLAIN) ../include/pgapack.h $(INSTALL_INCDIR)/
	$(INSTALL_PLAIN) ../include/pgapackf.h $(INSTALL_INCDIR)/
ifeq (${MPI},serial)
	$(INSTALL) -d $(INSTALL_INCDIR)/fakempi
	$(INSTALL_PLAIN) ../fakempi/mpi.h $(INSTALL_INCDIR)/fakempi/mpi.h
endif
	$(INSTALL) -m 0644 $(PCFILE) $(INSTALL_LIBDIR)/pkgconfig/
	ldconfig

# pkgconf config file
$(PCFILE):
	echo Name: PGAPack                                    > $@
	echo Description: Parallel Genetic Algorithm Package >> $@
	echo Version: $(VERSION)                             >> $@
	echo $(REQUIRES)                                     >> $@
	echo prefix=$(PREFIX)                                >> $@
	echo 'includedir=$${prefix}/include/pgapack-'$(MPI)  >> $@
	echo 'libdir=$${prefix}/'$(REL_LIBDIR)               >> $@
	echo 'Libs: -L$${libdir} -lpgapack-'$(MPI)           >> $@
	echo 'Cflags: -I$${includedir}' $(FAKEMPI_OPT)       >> $@

# Always re-create $(PCFILE)
.PHONY: clean $(PCFILE)

# vim: set noet :
