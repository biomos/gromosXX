INCLUDES = -I$(top_srcdir)/src

include_HEADERS =

bin_PROGRAMS =			md

md_SOURCES =			md.cc
md_LDADD = 			../src/libgroxx.la

EXTRA_DIST =			BUILD_NUMBER

autobuildnr :
	echo "#define BUILD_DATE "\""`date`"\" > BUILD_NUMBER.new
	awk '/BUILD_NUMBER/{ printf ("#define BUILD_NUMBER %d\n", $$3+1);}' ${srcdir}/BUILD_NUMBER >> BUILD_NUMBER.new
	if eval "test -w $(srcdir)/BUILD_NUMBER"; then \
		cp BUILD_NUMBER.new $(srcdir)/BUILD_NUMBER; \
	else \
		rm BUILD_NUMBER.new; \
	fi



