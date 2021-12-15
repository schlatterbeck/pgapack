all: classic dejong example maxbit maxchar maxint name namefull constr

%: %.f
	$(FC) $(FFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f classic dejong example maxbit maxchar maxint name namefull *.o

