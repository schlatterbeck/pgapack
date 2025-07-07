PROGS=testprog testprogde testprogdertr testprogmax
GCNO = $(PROGS:=.gcno)
GCDA = $(PROGS:=.gcda)
ROBJ = $(PROGS:=.o)

all: $(PROGS)

COMPILE.f = @$(FC) $(FFLAGS) -c

OBJS = grdfcn77.o hesfcn77.o objfcn77.o ssqfcn77.o ssqjac77.o vecfcn77.o \
       vecjac77.o

tprobs.a: $(OBJS)
	ar rv tprobs.a $(OBJS)
	$(RANLIB) tprobs.a

%.o: %.f
	$(FC) $(FFLAGS) -c $<

testprog: testprog.o tprobs.a
	$(FC) testprog.o tprobs.a $(LDFLAGS) -o testprog

testprogmax: testprogmax.o tprobs.a
	$(FC) testprogmax.o tprobs.a $(LDFLAGS) -o testprogmax

testprogde: testprogde.o tprobs.a
	$(CC) testprogde.o tprobs.a $(LDFLAGS) -o testprogde

testprogdertr: testprogdertr.o tprobs.a
	$(CC) testprogdertr.o tprobs.a $(LDFLAGS) -o testprogdertr

clean:
	rm -f $(OBJS) $(ROBJ) tprobs.a $(PROGS) $(GCDA) $(GCNO)

# vim: set noet :
