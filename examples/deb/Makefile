OBJS = deb00.o optimize.o

all: optimize

optimize: optimize.o $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o optimize

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f $(OBJS) optimize
