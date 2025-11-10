AM_CXXFLAGS = @MY_CXXFLAGS@
noinst_LTLIBRARIES = libfit.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = PositionUtils.h \
	Reference.h \
	RotationalFit.h \
	TranslationalFit.h \
	FastRotationalFit.h \
	AtomDistances.h

libfit_la_SOURCES = PositionUtils.cc \
	Reference.cc \
	RotationalFit.cc \
	TranslationalFit.cc \
	FastRotationalFit.cc \
	AtomDistances.cc

check_PROGRAMS = PositionUtils \
	Reference \
	RotationalFit \
	TranslationalFit

LDADD = ../libgromos.la

PositionUtils_SOURCES = PositionUtils.t.cc
Reference_SOURCES = Reference.t.cc
RotationalFit_SOURCES = RotationalFit.t.cc
TranslationalFit_SOURCES = TranslationalFit.t.cc


PositionUtils_LDADD =  \
	../fit/libfit.la \
	../utils/libutils.la \
	../bound/libbound.la \
	../gio/libgio.la \
	../gcore/libgcore.la \
	../args/libargs.la \
	../gmath/libgmath.la \
	$(GSL_LIB)
