## Template for automake for libgromos
AM_CXXFLAGS = @MY_CXXFLAGS@
SUBDIRS = gromos gmath gcore gio bound utils fit args pb .

EXTRA_DIST = Doxygen.h

lib_LTLIBRARIES = libgromos.la gromosplugin.la

libgromos_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0 $(GSL_LDFLAGS)
libgromos_la_SOURCES =
libgromos_la_LIBADD =	gmath/libgmath.la \
			gcore/libgcore.la \
			gio/libgio.la \
			bound/libbound.la \
			utils/libutils.la \
			fit/libfit.la \
			args/libargs.la \
			pb/libpb.la \
			$(GSL_LIB)

gromosplugin_la_LDFLAGS = -version-info $(GROMOS_MAJOR_VERSION):$(GROMOS_MICRO_VERSION):0 $(GSL_LDFLAGS) -module
gromosplugin_la_SOURCES =
gromosplugin_la_LIBADD =	gmath/libgmath.la \
			gcore/libgcore.la \
			gio/libgio.la \
			bound/libbound.la \
			utils/libutils.la \
			fit/libfit.la \
			args/libargs.la \
			pb/libpb.la \
			$(GSL_LIB)
