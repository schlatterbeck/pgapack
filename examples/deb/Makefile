OBJS = deb0.o deb1.o deb2.o deb3.o deb4.o deb5.o deb6.o deb7.o \
       deb8.o deb9.o optimize.o

all: optimize

optimize: optimize.o $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o optimize

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f $(OBJS) optimize

# vim: set noet :
