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
ifeq (${SHAREDLIBS},yes)
SHLIBFILE = lib$(PGA_LIB)-$(MPI).so.$(VERSION)
SHLIB = $(PGA_LIB_DIR)/$(SHLIBFILE)
SHARED = -shared -Wl,-soname,lib$(PGA_LIB)-$(MPI).so.$(MAJOR)
SHLIB_INSTALL = $(INSTALL) $(SHLIB) $(INSTALL_LIBDIR) && \
    ln -sf $(SHLIBFILE) $(INSTALL_LIBDIR)/lib$(PGA_LIB)-$(MPI).so
else
SHLIB_INSTALL = true
endif

all: $(PGA_LIB_DIR)/lib$(PGA_LIB).a $(SHLIB) $(PCFILE)

$(PGA_LIB_DIR)/lib$(PGA_LIB).a: $(OBJS)
	ar rv $(PGA_LIB_DIR)/lib$(PGA_LIB).a $(OBJS)
	$(RANLIB) $(PGA_LIB_DIR)/lib$(PGA_LIB).a

ifeq (${SHAREDLIBS},yes)
$(SHLIB): $(OBJS)
	$(CC) $(SHARED) $(OBJS) -o $(SHLIB)
endif

$(PGA_LIB_DIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(PCFILE)

install: all
	$(INSTALL) -d $(INSTALL_LIBDIR)/pkgconfig
	$(INSTALL) -d $(PREFIX)/include/pgapack-$(MPI)
	$(SHLIB_INSTALL)
	$(INSTALL) -T $(PGA_LIB_DIR)/lib$(PGA_LIB).a \
	    $(INSTALL_LIBDIR)/lib$(PGA_LIB)-$(MPI).a
	$(INSTALL) ../include/pgapack.h $(PREFIX)/include/pgapack-$(MPI)/
	$(INSTALL) ../include/pgapackf.h $(PREFIX)/include/pgapack-$(MPI)/
	$(INSTALL) -m 0644 $(PCFILE) $(INSTALL_LIBDIR)/pkgconfig/
	ldconfig

# pkgconf config file
$(PCFILE): Makefile
	echo Name: PGAPack                                    > $@
	echo Description: Parallel Genetic Algorithm Package >> $@
	echo Version: $(VERSION)                             >> $@
	echo Requires: mpi                                   >> $@
	echo prefix=$(PREFIX)                                >> $@
	echo 'includedir=$${prefix}/include/pgapack-'$(MPI)  >> $@
	echo 'libdir=$${prefix}/'$(REL_LIBDIR)               >> $@
	echo 'Libs: -L$${libdir} -lpgapack-'$(MPI)           >> $@
	echo 'Cflags: -I$${includedir}'                      >> $@

.PHONY: clean

# vim: set noet :
