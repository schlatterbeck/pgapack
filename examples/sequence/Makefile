OBJS = optimize.o

all: optimize magic

optimize: optimize.o
	$(CC) $(CFLAGS) optimize.o $(LDFLAGS) -o optimize

magic: magic.o
	$(CC) $(CFLAGS) magic.o $(LDFLAGS) -o magic

%.o: %.c optimize.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJS) optimize magic

# vim: set noet :
