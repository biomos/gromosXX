SUBDIRS = src program

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

EXTRA_DIST = VERSION \
	doc/doxygen.conf.in \
	doc/newinput.txt \
	doc/arguments.doxy \
	doc/arguments.html \
	doc/input.doxy \
	doc/input.html \
	doc/install.doxy \
	doc/install.html \
	doc/main.doxy \
	doc/output.doxy \
	doc/output.html \
	doc/posres.doxy \
	doc/posres.html \
	doc/pttopo.doxy \
	doc/pttopo.html \
	doc/release.doxy \
	doc/release.html \
	doc/tail.html \
	packaging/solaris.depend \
	packaging/debian.control \
	packaging/copyright

