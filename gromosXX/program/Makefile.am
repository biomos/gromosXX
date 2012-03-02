AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src

include_HEADERS =

bin_PROGRAMS = md md_mpi repex_mpi eds_2box eds_2box_mpi

md_SOURCES = md.cc
md_mpi_SOURCES = md_mpi.cc
repex_mpi_SOURCES = repex_mpi.cc
eds_2box_SOURCES = eds_2box.cc
eds_2box_mpi_SOURCES = eds_2box_mpi.cc

LDADD = ../src/libmdpp.la $(extra_lib)

EXTRA_DIST = BUILD_NUMBER

autobuildnr :
	echo "#define BUILD_DATE "\""`date`"\" > BUILD_NUMBER.new
	awk '/BUILD_NUMBER/{ printf ("#define BUILD_NUMBER %d\n", $$3+1);}' ${srcdir}/BUILD_NUMBER >> BUILD_NUMBER.new
	if eval "test -w $(srcdir)/BUILD_NUMBER"; then \
		cp BUILD_NUMBER.new $(srcdir)/BUILD_NUMBER; \
	else \
		rm BUILD_NUMBER.new; \
	fi



