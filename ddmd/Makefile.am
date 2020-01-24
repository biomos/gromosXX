SUBDIRS = .

docdir = $(prefix)/doc/$(PACKAGE)-$(VERSION)

doc:	doc/doxygen.conf.in
	sed "s/<SRCDIR>/$(top_srcdir)/g" $(top_srcdir)/doc/doxygen.conf.in > doxygen.conf
	mkdir -p doc
	doxygen doxygen.conf

doc-install:	doc
	install -d $(docdir)
	cp -r doc/html $(docdir)


bin_PROGRAMS =	ddmd

ddmd_SOURCES =  md/ddmd.cc \
		simulation/multibath.cc \
		io/blockinput.cc \
		io/instream.cc \
		io/message.cc \
		io/usage.cc \
		io/argument.cc \
		io/read_input.cc \
		io/in_parameter.cc \
		simulation/simulation.cc \
		gmath.cc
