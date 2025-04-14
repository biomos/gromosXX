AM_CXXFLAGS = @MY_CXXFLAGS@
ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects
SUBDIRS = src program contrib data

docdir = $(prefix)/doc/$(PACKAGE)-$(VERSION)

build:
	cd program ; make autobuildnr; cd ..
	make

doc:	doc/doxygen.conf.in
	sed "s/<SRCDIR>/$(top_srcdir)/g" $(top_srcdir)/doc/doxygen.conf.in > doxygen.conf
	mkdir -p  $(top_srcdir)/doc/snippets
	python $(top_srcdir)/doc/generate_snippets.py $(top_srcdir)/src/io/parameter/in_parameter.cc $(top_srcdir)/src/io/topology/in_*.cc > $(top_srcdir)/doc/snippets/snippets.cc
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
	doc/main.doxy\
    data/ene_ana.md++.lib\
	contrib/schnet_v1.py
	contrib/schnet_v2.py
