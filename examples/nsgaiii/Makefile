OBJS = optimize.o dtlz1.o dtlz2.o dtlz3.o dtlz4.o scaled_dtlz1.o \
    scaled_dtlz2.o convex_dtlz2.o neg_dtlz2.o c1_dtlz1.o c1_dtlz3.o \
    c2_dtlz2.o c2_convex_dtlz2.o

all: optimize

optimize: optimize.o $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o optimize

%.o: %.c constraint.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f $(OBJS) optimize
