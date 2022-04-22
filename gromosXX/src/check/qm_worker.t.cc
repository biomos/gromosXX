#include <algorithm>

#include "qm_worker.t.h"

#include "../io/argument.h"
#include "../io/read_input.h"

#include "../util/usage.h"


int QM_Worker_Test::init_simulation() {
  // supported files
  std::vector<std::string> qualifiers = {"topo", "conf", "input", "qmmm", "trc", "tre", "fin"};
  util::Known knowns;
  for (std::string& e : qualifiers) { knowns << e; }
  // process arguments
  std::for_each(qualifiers.begin(), qualifiers.end(), [](std::string& s) {s.insert(0, "@"); });
  io::Argument args;
  char* c_style_args[] = {
    &this->binary_name[0],
    &qualifiers[0][0],
    &this->topology_file[0],
    &qualifiers[1][0],
    &this->configuration_file[0],
    &qualifiers[2][0],
    &this->simulation_file[0],
    &qualifiers[3][0],
    &this->qmmm_file[0],
    &qualifiers[4][0],
    &this->coordinate_trajectory_file[0],
    &qualifiers[5][0],
    &this->energy_trajectory_file[0],
    &qualifiers[6][0],
    &this->final_configuration_file[0],
  };
  args.parse(15, c_style_args, knowns);

  // initialize the simulation
  io::read_input(args, this->topo, this->conf, this->sim, this->md);
  this->md.init(this->topo, this->conf, this->sim);

  return EXIT_SUCCESS;
}

QM_Worker_Test::~QM_Worker_Test() {

}