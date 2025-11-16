AM_CXXFLAGS = @MY_CXXFLAGS@
noinst_LTLIBRARIES = libgio.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS =  Ginstream.h\
	gzstream.h\
	InTopology.h\
	OutTopology.h\
	OutGromacs.h\
	InG96.h\
	OutCoordinates.h\
	OutG96.h\
	OutG96S.h\
	OutPdb.h\
	Outvmdam.h\
	OutCif.h\
	InParameter.h\
	InBuildingBlock.h\
        InLinkingBlock.h\
	InPtTopology.h\
	OutPtTopology.h\
	StringTokenizer.h\
	InCIF.h\
	InIACElementNameMapping.h\
	InBFactorOccupancy.h\
	OutBuildingBlock.h\
	InPDB.h\
	InAmberTopology.h\
	InChargeGroups.h

libgio_la_SOURCES = Ginstream.cc\
	gzstream.cc\
	InG96.cc\
	InTopology.cc\
	InParameter.cc\
	InBuildingBlock.cc\
        InLinkingBlock.cc\
	InPtTopology.cc\
	OutTopology.cc\
	OutPtTopology.cc\
	OutGromacs.cc\
	OutCoordinates.cc\
	OutG96.cc\
	OutG96S.cc\
	OutPdb.cc\
	Outvmdam.cc\
	OutCif.cc\
	StringTokenizer.cc\
	InCIF.cc\
	InIACElementNameMapping.cc\
	InBFactorOccupancy.cc\
	OutBuildingBlock.cc\
    	InPDB.cc\
    	InAmberTopology.cc\
	InChargeGroups.cc

check_PROGRAMS = InG96\
	InTopology\
	InParameter\
	InBuildingBlock\
	OutGromacs\
	OutG96S\
	OutPdb\
	Outvmdam

AM_LDFLAGS = $(GSL_LDFLAGS)
LDADD = libgio.la \
	../gcore/libgcore.la \
	../gmath/libgmath.la \
	../args/libargs.la \
	../utils/libutils.la \
	$(GSL_LIB)

InG96_SOURCES = InG96.t.cc
InTopology_SOURCES = InTopology.t.cc
InParameter_SOURCES = InParameter.t.cc
InBuildingBlock_SOURCES = InBuildingBlock.t.cc
OutGromacs_SOURCES = OutGromacs.t.cc
OutG96S_SOURCES = OutG96S.t.cc
OutPdb_SOURCES = OutPdb.t.cc
Outvmdam_SOURCES = Outvmdam.t.cc
