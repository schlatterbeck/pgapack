ALL=classic dejong example maxbit maxchar maxint name namefull constr
all: ${ALL}

%: %.f
	$(FC) $(FFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f ${ALL} *.o
