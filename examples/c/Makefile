PROG=classic dejong example maxbit maxchar maxint name namefull udtstr \
    leastsquare deb10 maxint-de maxbit-hc
GCNO=$(PROG:=.gcno)
GCDA=$(PROG:=.gcda)

all: $(PROG)

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f $(PROG) $(GCNO) $(GCDA) *.o

# vim: set noet :
