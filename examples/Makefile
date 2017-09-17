all:
	@echo "Making C examples"
	@make C
	@echo "Making fortran examples"
	@make Fortran

C:
	@cd c; $(MAKE)

Fortran:
	@if test -n "$(FC)" ; then \
		cd fortran && $(MAKE); cd ..; \
		cd mgh && $(MAKE); \
	else \
		echo "Fortran compiler unavailable"; \
	fi

clean:
	@cd c; $(MAKE) clean
	@cd fortran; $(MAKE) clean
	@cd mgh; $(MAKE) clean

