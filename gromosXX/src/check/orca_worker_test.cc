/**
 * @file orca_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"
#include "../algorithm/algorithm.h"
#include "../algorithm/algorithm/algorithm_sequence.h"

#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../io/configuration/out_configuration.h"

#include "check.h"
#include "test_parameters.h"

#include "qm_worker_test.h"
#include "orca_worker_test.h"
#include "orca_worker_test_electrostatic.h"
#include "orca_worker_test_mechanical.h"

#include "../math/gmath.h"

#include "../interaction/qmmm/orca_worker.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

namespace testing {

void Orca_Worker_Test::init_results_binary_name() {
  // binary and worker name
  results_.strs_[Key::binary] = "/home/fpultar/opt/orca-5.0.3/orca"; // modify this
  results_.strs_[Key::qm_worker_name] = "Orca Worker";
}

void Orca_Worker_Test::init_results_units() {
  // units and conversion factors
  results_.doubles_[Key::unit_factor_length] = 0.1;
  results_.doubles_[Key::unit_factor_energy] = 2625.5;
  results_.doubles_[Key::unit_factor_force] = 49641;
  results_.doubles_[Key::unit_factor_charge] = 1.0;
}

void Orca_Worker_Test::init_results_elements() {
  // element symbol and iac atom mapping
  results_.strs_[Key::symbol_h] = "H";
  results_.strs_[Key::symbol_c] = "C";
  results_.strs_[Key::symbol_o] = "O";
  results_.strs_[Key::symbol_p] = "P";
  results_.strs_[Key::symbol_pd] = "Pd";

  results_.ints_[Key::element_1] = 1;
  results_.ints_[Key::element_6] = 6;
  results_.ints_[Key::element_8] = 8;
  results_.ints_[Key::element_15] = 15;
  results_.ints_[Key::element_20] = 20;
  results_.ints_[Key::element_26] = 26;
  results_.ints_[Key::element_29] = 29;
  results_.ints_[Key::element_46] = 46;
}

} // namespace testing