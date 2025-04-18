PROG=classic dejong example maxbit maxchar maxint name namefull udtstr \
    leastsquare deb10 maxint-de maxbit-hc

all: $(PROG)

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f $(PROG) *.o

# vim: set noet :
