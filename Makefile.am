AM_CXXFLAGS = @MY_CXXFLAGS@
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src program contrib

docdir = $(prefix)/doc/$(PACKAGE)-$(VERSION)

build:
	cd program ; make autobuildnr; cd ..
	make

doc:	doc/doxygen.conf.in
	sed "s/<SRCDIR>/$(top_srcdir)/g" $(top_srcdir)/doc/doxygen.conf.in > doxygen.conf
	mkdir -p doc
	doxygen doxygen.conf

doc-install:	doc
	install -d $(docdir)
	cp -r doc/html $(docdir)

EXTRA_DIST = VERSION\
	contrib/tabulate_spc.cc\
	doc/doxygen.conf.in\
	doc/install.doxy\
	doc/install.html\
	doc/main.doxy
