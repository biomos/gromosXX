AM_CXXFLAGS = @MY_CXXFLAGS@
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
	contrib/repana.cc\
	contrib/rep_rewrite.cc\
	contrib/tabulate_spc.cc\
	doc/doxygen.conf.in\
	doc/newinput.txt\
	doc/arguments.doxy\
	doc/arguments.html\
	doc/input.doxy\
	doc/input.html\
	doc/install.doxy\
	doc/install.html\
	doc/main.doxy\
	doc/output.doxy\
	doc/output.html\
	doc/posres.doxy\
	doc/posres.html\
	doc/release.doxy\
	doc/release.html\
	doc/tail.html\
	packaging/solaris.depend\
	packaging/debian.control \
	packaging/debian.depend\
	packaging/copyright\
	nightly/nightly.sh\
	nightly/prepare.sh\
	nightly/distcheck.sh\
	nightly/options.gromos++\
	nightly/options.gromosXX\
	nightly/options.all\
	nightly/options.SunOS_i86pc\
	nightly/options.Linux_i686\
	nightly/options.SunOS_sun4u\
	nightly/pkg_solaris.sh\
	nightly/pkg_solaris_fakeroot.sh\
	nightly/pkg_debian.sh\
	nightly/pkg_debian_fakeroot.sh


