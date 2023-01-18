# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SOURCE        := ${wildcard ../../source/*.c}
INCLUDE       := ../../include/pgapack.h
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     := .
BUILDDIR      := _build
MAN_PATH1     := $(BUILDDIR)/man/man3a
MAN_PATH2     := man/man3
MAN_STAGE1     = ${wildcard $(MAN_PATH1)/*.3a}
MAN_STAGE1_3   = ${patsubst %.3a,%.3,$(MAN_STAGE1)}
MAN_STAGE2     = ${subst $(MAN_PATH1),$(MAN_PATH2),$(MAN_STAGE1_3)}

all: html fixedman latexpdf
	@echo Done

html: .stamp-doxygen

man: .stamp-doxygen

latex: .stamp-doxygen

latexpdf: .stamp-doxygen

fixedman: man man/man3
	$(MAKE) fixedman2

fixedman2: $(MAN_STAGE2)
	cp $(MAN_PATH2)/_PGAGetEvaluation.3 $(MAN_PATH2)/PGAGetEvaluation.3
	cp $(MAN_PATH2)/_PGASetEvaluation.3 $(MAN_PATH2)/PGASetEvaluation.3
	@echo Stage2 done

.stamp-doxygen: $(SOURCE) $(INCLUDE) Doxyfile
	doxygen
	touch .stamp-doxygen

Doxyfile:
	@echo Doxyfile

$(MAN_PATH2):
	mkdir -p $(MAN_PATH2)

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile clean fixedman2 fixedman $(MAN_STAGE1)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf _doxygen api _static  _templates man.py man/man3
	rm -f .stamp-*
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

$(MAN_PATH2)/%.3: $(MAN_PATH1)/%.3a
	fixman < $< > $@

# vim: set noet :
