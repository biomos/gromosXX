AM_CXXFLAGS = @MY_CXXFLAGS@
noinst_LTLIBRARIES = libbound.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS     = Boundary.h \
		  Vacuum.h \
		  TruncOct.h \
		  RectBox.h \
                  Triclinic.h

libbound_la_SOURCES = Boundary.cc\
                     Vacuum.cc\
                     TruncOct.cc\
                     RectBox.cc \
                     Triclinic.cc

check_PROGRAMS = TruncOct RectBox

LDADD = ../libgromos.la

TruncOct_SOURCES = TruncOct.t.cc 
RectBox_SOURCES = RectBox.t.cc
