OBJS = optimize.o dtlz1.o dtlz2.o

all: optimize

optimize: optimize.o $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o optimize

%.o: %.c constraint.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f $(OBJS) optimize
