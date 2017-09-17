COMPILE.f = @$(FC) $(FFLAGS) -c

OBJS = grdfcn77.o hesfcn77.o objfcn77.o ssqfcn77.o ssqjac77.o vecfcn77.o \
       vecjac77.o

tprobs.a: $(OBJS)
	ar rv tprobs.a $(OBJS)
	$(RANLIB) tprobs.a

%.o: %.f
	$(FC) $(FFLAGS) -c $<

testprog: testprog.f tprobs.a
	$(FC) $(FFLAGS) testprog.f tprobs.a $(LDFLAGS)

clean:
	rm -f $(OBJS) tprobs.a testprog
