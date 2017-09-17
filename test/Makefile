all: instverf 

%: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean: 
	rm -f instverf instverf.o

.PHONY: all clean
