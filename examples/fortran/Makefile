ALL=classic dejong example maxbit maxchar maxint name namefull constr
GCNO = $(ALL:=.gcno)
GCDA = $(ALL:=.gcda)

all: ${ALL}

%: %.f
	$(FC) $(FFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f $(ALL) $(GCNO) $(GCDA) *.o

# vim: set noet :
