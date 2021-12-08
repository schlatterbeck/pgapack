all: C Fortran

C:
	$(MAKE) -C c
	$(MAKE) -C deb

ifeq (,${FC})

Fortran:
	@echo Fortran compiler unavailable

else

Fortran:
	$(MAKE) -C fortran
	$(MAKE) -C mgh

endif

clean:
	$(MAKE) -C c       clean
	$(MAKE) -C fortran clean
	$(MAKE) -C mgh     clean
	$(MAKE) -C deb     clean

.PHONY: all C Fortran clean
