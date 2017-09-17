
all: classic dejong example maxbit maxchar maxint name namefull udtstr

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f classic dejong example maxbit maxchar maxint name namefull \
	       udtstr *.o
