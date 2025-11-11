AM_CXXFLAGS = @MY_CXXFLAGS@
noinst_LTLIBRARIES = libutils.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = parse.h \
	Rmsd.h \
	Neighbours.h \
	VirtualAtom.h \
	Noe.h \
	AtomSpecifier.h \
	AtomicRadii.h \
	CremerPople.h \
	VectorSpecifier.h \
	PropertyContainer.h \
	Property.h \
	Energy.h \
	TrajArray.h \
	CheckTopo.h \
	Hbond_calc.h \
	Hbond_calc_2c.h \
	Hbond_calc_3c.h \
	Hbond_calc_bridges.h \
	Hbond.h \
	Dssp.h \
	EnergyTraj.h \
	RestrTraj.h \
	DipTraj.h \
	SimplePairlist.h \
	SoluteWeightedDistance.h \
	FfExpert.h \
	FfExpertGraph.h \
	make_top.h \
	Temperature.h \
	Value.h \
	ExpressionParser.h \
	ExpressionParser.cc \
	groTime.h \
	RdcFuncs.h \
	JvalFuncs.h \
	CommandLine.h \
	debug.h\
	AminoAcid.h\
	RDF.h\
	NeutronScattering.h\
    CubeSystem.hcc\
	Disicl.h\
	Gch.h\
	IntegerInputParser.h\
	StringOps.h

libutils_la_SOURCES = parse.cc \
	Rmsd.cc \
	Neighbours.cc \
	VirtualAtom.cc \
	Noe.cc \
	AtomSpecifier.cc \
	AtomSpecifier_parse.cc \
	AtomicRadii.cc \
	CremerPople.cc \
	VectorSpecifier.cc \
	PropertyContainer.cc \
	Property.cc \
	Energy.cc \
	TrajArray.cc \
	CheckTopo.cc \
	Hbond_calc.cc \
	Hbond_calc_2c.cc \
	Hbond_calc_3c.cc \
	Hbond_calc_bridges.cc \
	Hbond.cc \
	Dssp.cc \
	EnergyTraj.cc \
	RestrTraj.cc \
	DipTraj.cc \
	SimplePairlist.cc \
	SoluteWeightedDistance.cc \
	FfExpert.cc \
	FfExpertGraph.cc \
	groTime.cc \
	RdcFuncs.cc \
	Temperature.cc \
	JvalFuncs.cc \
	CommandLine.cc \
	debug.cc\
	AminoAcid.cc\
	RDF.cc\
	NeutronScattering.cc\
	Disicl.cc\
	Gch.cc\
	IntegerInputParser.cc\
	StringOps.cc

check_PROGRAMS = ExpressionParser \
	Rmsd \
	Neighbours \
	VirtualAtom \
	Noe \
	AtomSpecifier \
	PropertyContainer \
	Energy \
	CheckTopo \
	SimplePairlist \
	FfExpert


LDADD = ../libgromos.la

ExpressionParser_SOURCES = ExpressionParser.t.cc
Rmsd_SOURCES = Rmsd.t.cc
Neighbours_SOURCES = Neighbours.t.cc
VirtualAtom_SOURCES = VirtualAtom.t.cc
Noe_SOURCES = Noe.t.cc
AtomSpecifier_SOURCES = AtomSpecifier.t.cc
PropertyContainer_SOURCES = PropertyContainer.t.cc
Energy_SOURCES = Energy.t.cc
CheckTopo_SOURCES = CheckTopo.t.cc
SimplePairlist_SOURCES = SimplePairlist.t.cc
FfExpert_SOURCES = FfExpert.t.cc

AM_LDFLAGS = $(GSL_LDFLAGS)

VirtualAtom_LDADD = libutils.la \
	../bound/libbound.la \
	../fit/libfit.la \
	../gio/libgio.la \
	../gcore/libgcore.la \
	../args/libargs.la \
	../gmath/libgmath.la \
	$(GSL_LIB)

FfExpert_LDADD = libutils.la \
	../bound/libbound.la \
	../fit/libfit.la \
	../gio/libgio.la \
	../gcore/libgcore.la \
	../args/libargs.la \
	../gmath/libgmath.la \
	$(GSL_LIB)

ExpressionParser_LDADD = libutils.la \
	../gio/libgio.la \
	../gcore/libgcore.la \
	../args/libargs.la \
	../gmath/libgmath.la \
	../bound/libbound.la \
	../fit/libfit.la \
	$(GSL_LIB)

