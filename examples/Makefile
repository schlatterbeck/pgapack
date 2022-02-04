all: C Fortran

C:
	$(MAKE) -C c
	$(MAKE) -C deb
	$(MAKE) -C nsgaii
	$(MAKE) -C nsgaiii
	$(MAKE) -C sequence

ifeq (,${FC})

Fortran:
	@echo Fortran compiler unavailable

else

Fortran:
	$(MAKE) -C fortran
	$(MAKE) -C mgh

endif

clean:
	$(MAKE) -C c        clean
	$(MAKE) -C fortran  clean
	$(MAKE) -C mgh      clean
	$(MAKE) -C deb      clean
	$(MAKE) -C nsgaii   clean
	$(MAKE) -C nsgaiii  clean
	$(MAKE) -C sequence clean

.PHONY: all C Fortran clean
