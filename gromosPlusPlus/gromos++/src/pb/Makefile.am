AM_CXXFLAGS = @MY_CXXFLAGS@
noinst_LTLIBRARIES = libpb.la

gincludedir = $(includedir)/gromos++

ginclude_HEADERS = FDPoissonBoltzmann.h\
   		   FDPoissonBoltzmann_ICCG_NPBC.h\
   		   FDPoissonBoltzmann_ICCG_PBC.h\
   		   PB_Parameters.h\
		   FFTBoundaryCondition.h\
		   FFTChargeDipole.h\
	           FFTChargeDipole_LS.h\
                   FFTChargeDipole_RF.h\
                   FFTChargeShapingFunction.h\
                   FFTChargeShapingFunction_Hat.h\
                   FFTChargeShapingFunction_Parabola.h\
                   FFTDipoleDipole.h\
                   FFTDipoleDipole_LS.h\
                   FFTDipoleDipole_RF.h\
                   FFTGridType.h\
                   FFTInsideOutside.h\
                   FFTInteractionTypeCodes.h\
                   FFTPoisson.h\
                   FFTPoissonIterator.h\
                   FFTVacuumField.h\
                   FFTVacuumField_LS.h\
                   FFTVacuumField_RF.h\
                   Ewald_edir.h

libpb_la_SOURCES = FDPoissonBoltzmann.cc\
   		   FDPoissonBoltzmann_ICCG_NPBC.cc\
   		   FDPoissonBoltzmann_ICCG_PBC.cc\
   		   PB_Parameters.cc\
		   FFTBoundaryCondition.cc\
		   FFTChargeDipole.cc\
	           FFTChargeDipole_LS.cc\
                   FFTChargeDipole_RF.cc\
                   FFTChargeShapingFunction.cc\
                   FFTChargeShapingFunction_Hat.cc\
                   FFTChargeShapingFunction_Parabola.cc\
                   FFTDipoleDipole.cc\
                   FFTDipoleDipole_LS.cc\
                   FFTDipoleDipole_RF.cc\
                   FFTGridType.cc\
                   FFTInsideOutside.cc\
                   FFTInteractionTypeCodes.cc\
                   FFTPoisson.cc\
                   FFTPoissonIterator.cc\
                   FFTVacuumField.cc\
                   FFTVacuumField_LS.cc\
                   FFTVacuumField_RF.cc\
                   Ewald_edir.cc


LDADD = libpb.la

