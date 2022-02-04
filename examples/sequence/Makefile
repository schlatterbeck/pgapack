OBJS = optimize.o

all: optimize

optimize: optimize.o $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o optimize

%.o: %.c optimize.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f $(OBJS) optimize
