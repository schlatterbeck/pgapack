OBJS = optimize.o magic.o
GCNO = $(OBJS:.o=.gcno)
GCDA = $(OBJS:.o=.gcda)

all: optimize magic

optimize: optimize.o
	$(CC) $(CFLAGS) optimize.o $(LDFLAGS) -o optimize

magic: magic.o
	$(CC) $(CFLAGS) magic.o $(LDFLAGS) -o magic

%.o: %.c optimize.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJS) $(GCNO) $(GCDA) optimize magic

# vim: set noet :
