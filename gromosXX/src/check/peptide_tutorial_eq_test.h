/**
 * @file peptide_tutorial_eq_test.h
 * Class to perform tests on the equilibration step of the peptide tutorial
 */
#ifndef INCLUDED_PEPTIDE_TUTORIAL_EQ_TEST_H
#define	INCLUDED_PEPTIDE_TUTORIAL_EQ_TEST_H

#include "../stdheader.h"
#include "simulation_test.h"

namespace testing {
  
class Peptide_Tutorial_Eq_Test : public Simulation_Test {

public:

  Peptide_Tutorial_Eq_Test(const std::string& name) : Simulation_Test(Parameter(name, "tutorial_test", "src/check/data/tutorial/eq")) {}

};
  
}

#endif /* INCLUDED_PEPTIDE_TUTORIAL_EQ_TEST_H */