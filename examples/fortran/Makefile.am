check_PROGRAMS= classic dejong example maxbit maxchar name namefull
AM_FFLAGS=-I../../include -DWL=32

classic_SOURCES=classic.f
dejong_SOURCES=dejong.f
example_SOURCES=example.f
maxbit_SOURCES=maxbit.f
maxchar_SOURCES=maxchar.f
name_SOURCES=name.f
namefull_SOURCES=namefull.f

classic_LDADD= -L../../source -lpgapack -lmpi -lm $(FLIBS)
dejong_LDADD= -L../../source -lpgapack -lmpi -lm $(FLIBS)
example_LDADD= -L../../source -lpgapack -lmpi -lm $(FLIBS)
maxbit_LDADD= -L../../source -lpgapack -lmpi -lm $(FLIBS)
maxchar_LDADD= -L../../source -lpgapack -lmpi -lm $(FLIBS)
name_LDADD= -L../../source -lpgapack -lmpi -lm $(FLIBS)
namefull_LDADD= -L../../source -lpgapack -lmpi -lm $(FLIBS)
