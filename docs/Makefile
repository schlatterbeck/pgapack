SOURCE=user_guide.tex user_guide.bib tpage.sty
EPS=generation.eps refpoints.eps

PNGTOPPM=pngtopnm -mix -background white
PNMTOPS=pnmtops
PNMTOPS_OPT=-dpi 300 -equalpixels -noturn -rle

all: user_guide.pdf

user_guide.dvi: $(SOURCE) $(EPS) userguide_date.tex
	latex user_guide
	bibtex user_guide
	latex user_guide
	latex user_guide

userguide_date.tex: $(SOURCE) $(EPS)
	echo '\\def\\updated{%' > $@
	git log -n1 $(SOURCE) $(EPS) | grep '^Date:' | cut -d' ' -f5,8 >> $@
	echo '}' >> $@

user_guide.ps: user_guide.dvi
	dvips -t letter user_guide

user_guide.pdf: user_guide.ps
	ps2pdf user_guide.ps
	exiftool \
	-Title="Users Guide to the PGAPack Parallel Genetic Algorithm Library" \
	-Author="David Levine and contributors" \
	user_guide.pdf

%.ppm: %.png
	$(PNGTOPPM) $< > $@

%.eps: %.ppm
	$(PNMTOPS) $(PNMTOPS_OPT) $< > $@

clean:
	rm -f user_guide.aux user_guide.bbl user_guide.blg user_guide.dvi \
	    user_guide.idx user_guide.log user_guide.pdf user_guide.ps    \
	    user_guide.toc refpoints.eps refpoints.ppm userguide_date.tex \
	    user_guide.pdf_original


.PRECIOUS: user_guide.ps

# vim: set noet :
