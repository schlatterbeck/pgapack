SOURCE=user_guide.tex user_guide.bib tpage.sty

all: user_guide.pdf

user_guide.dvi: $(SOURCE)
	latex user_guide
	bibtex user_guide
	latex user_guide
	latex user_guide

user_guide.ps: user_guide.dvi
	dvips -t letter user_guide

user_guide.pdf: user_guide.ps
	ps2pdf user_guide.ps
	exiftool \
	-Title="Users Guide to the PGAPack Parallel Genetic Algorithm Library" \
	-Author="David Levine and contributors" \
        user_guide.pdf


clean:
	rm -f user_guide.aux user_guide.bbl user_guide.blg user_guide.dvi \
	    user_guide.idx user_guide.log user_guide.pdf user_guide.ps    \
	    user_guide.toc


.PRECIOUS: user_guide.ps
