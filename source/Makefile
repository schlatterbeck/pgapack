OBJS = $(PGA_LIB_DIR)/binary.o        \
      $(PGA_LIB_DIR)/char.o          \
      $(PGA_LIB_DIR)/cmdline.o       \
      $(PGA_LIB_DIR)/create.o        \
      $(PGA_LIB_DIR)/cross.o         \
      $(PGA_LIB_DIR)/debug.o         \
      $(PGA_LIB_DIR)/duplcate.o      \
      $(PGA_LIB_DIR)/evaluate.o      \
      $(PGA_LIB_DIR)/fitness.o       \
      $(PGA_LIB_DIR)/hamming.o       \
      $(PGA_LIB_DIR)/heap.o          \
      $(PGA_LIB_DIR)/integer.o       \
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

########################################################
#                    The Target                        #
########################################################

MAJOR = 1
MINOR = 0
ifeq (${SHAREDLIBS},yes)
SHLIB = $(PGA_LIB_DIR)/lib$(PGA_LIB)-$(MPI).so.$(MAJOR).$(MINOR)
SHARED = -shared -Wl,-soname,$(PGA_LIB)-$(MPI).so.$(MAJOR)
endif

all: $(PGA_LIB_DIR)/lib$(PGA_LIB).a $(SHLIB)

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
	rm -f $(OBJS)

.PHONY: clean
