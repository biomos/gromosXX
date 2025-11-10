# This is the File that generates Makefile.in
# through automake.
AM_CXXFLAGS = @MY_CXXFLAGS@
ACLOCAL_AMFLAGS = -I m4
docdir = $(prefix)/doc/$(PACKAGE)-$(VERSION)

SUBDIRS = src \
	programs \
	contrib \
	data

EXTRA_DIST =	VERSION \
	doc/doxygen.conf.in 

doc:	doc/doxygen.conf.in programs/0index.doxy contrib/0index_contrib.doxy
	sed "s/<SRCDIR>/${top_srcdir}/g" ${top_srcdir}/doc/doxygen.conf.in > doxygen.conf
	mkdir -p doc
	doxygen doxygen.conf

doc-install:	doc
	install -d $(docdir)
	cp -r doc/html $(docdir)
