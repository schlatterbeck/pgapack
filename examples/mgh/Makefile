PROGS=testprog testprogde testprogdertr testprogmax
all: $(PROGS)

COMPILE.f = @$(FC) $(FFLAGS) -c

OBJS = grdfcn77.o hesfcn77.o objfcn77.o ssqfcn77.o ssqjac77.o vecfcn77.o \
       vecjac77.o

tprobs.a: $(OBJS)
	ar rv tprobs.a $(OBJS)
	$(RANLIB) tprobs.a

%.o: %.f
	$(FC) $(FFLAGS) -c $<

testprog: testprog.f tprobs.a
	$(FC) $(FFLAGS) testprog.f tprobs.a $(LDFLAGS) -o testprog

testprogmax: testprogmax.f tprobs.a
	$(FC) $(FFLAGS) testprogmax.f tprobs.a $(LDFLAGS) -o testprogmax

testprogde: testprogde.c tprobs.a
	$(CC) $(CFLAGS) testprogde.c tprobs.a $(LDFLAGS) -o testprogde

testprogdertr: testprogdertr.c tprobs.a
	$(CC) $(CFLAGS) testprogdertr.c tprobs.a $(LDFLAGS) -o testprogdertr

clean:
	rm -f $(OBJS) tprobs.a $(PROGS)
