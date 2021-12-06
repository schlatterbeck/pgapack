PROG=classic dejong example maxbit maxchar maxint name namefull udtstr \
    leastsquare

all: $(PROG)

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f $(PROG) *.o
