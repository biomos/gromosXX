AM_CXXFLAGS = @MY_CXXFLAGS@
noinst_LTLIBRARIES = libgcore.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = Angle.h \
	AngleType.h \
	AtomPair.h \
	AtomTopology.h \
	MassType.h \
	Bond.h \
	BondType.h \
	Box.h \
	Remd.h \
	Constraint.h \
	CrossDihedral.h \
	Dihedral.h \
	DihedralType.h \
	Exclusion.h \
	GromosForceField.h \
	Improper.h \
	ImproperType.h \
	LJType.h \
	LJException.h \
	LJExceptionType.h \
	CGType.h \
	MoleculeTopology.h \
	Molecule.h \
	SolventTopology.h \
	Solvent.h \
	BuildingBlock.h \
	BbSolute.h \
        BbLink.h \
	System.h \
	LinearTopology.h \
	PtTopology.h \
        VirtualAtomType.h \
        VirtualAtoms.h \
	Weight.h

libgcore_la_SOURCES = Angle.cc \
	AngleType.cc \
	AtomPair.cc \
	AtomTopology.cc \
	MassType.cc \
	Bond.cc \
	BondType.cc \
	Box.cc \
	Remd.cc \
	Constraint.cc \
	CrossDihedral.cc \
	Dihedral.cc \
	DihedralType.cc \
	Exclusion.cc \
	GromosForceField.cc \
	Improper.cc \
	ImproperType.cc \
	LJType.cc \
	LJException.cc \
	LJExceptionType.cc \
	MoleculeTopology.cc \
	Molecule.cc \
	SolventTopology.cc \
	Solvent.cc \
	BuildingBlock.cc \
	BbSolute.cc \
        BbLink.cc \
	System.cc \
	LinearTopology.cc \
	PtTopology.cc \
        VirtualAtoms.cc \
        VirtualAtomType.cc \
	Weight.cc

check_PROGRAMS = AtomTopology \
	Exclusion \
	MoleculeTopology \
	Molecule \
	SolventTopology \
	Solvent \
	BbSolute

LDADD = libgcore.la \
	../args/libargs.la

AtomTopology_SOURCES = AtomTopology.t.cc
Exclusion_SOURCES = Exclusion.t.cc
MoleculeTopology_SOURCES = MoleculeTopology.t.cc
Molecule_SOURCES = Molecule.t.cc
SolventTopology_SOURCES = SolventTopology.t.cc
Solvent_SOURCES = Solvent.t.cc
BbSolute_SOURCES = BbSolute.t.cc
