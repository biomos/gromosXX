#include "test_parameters.h"

namespace testing {

Parameter orca_parameter_mechanical("Orca QM/MM Tests Mechanical Embedding", "orca_worker_test", "src/check/data/orca");

Results orca_results_mechanical;

Parameter orca_parameter_electrostatic("Orca QM/MM Tests Electrostatic Embedding", "orca_worker_test", "src/check/data/orca");

Results orca_results_electrostatic;

}