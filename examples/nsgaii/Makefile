OBJS = sch.o fon.o kur.o pol.o zdt1.o zdt2.o zdt3.o zdt4.o zdt6.o \
    constr.o srn.o tnk.o water.o rotated.o optimize.o

all: optimize

optimize: optimize.o $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o optimize

%.o: %.c constraint.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f $(OBJS) optimize
