EXTRA_DIST = Config.sh VERSION cvs2cl.pl doc/doxygen.conf doc/newinput.txt
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