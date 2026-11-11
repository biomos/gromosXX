#include "../../../stdheader.h"
#include "../../../algorithm/algorithm.h"
#include "../../../configuration/configuration.h"
#include "colvar_bias.h"
#include "colvar_state.h"

bool interaction::check_colvar_average_count(const configuration::Configuration &conf,
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs, bool perturbed)
{
  const auto &state = perturbed ? conf.special().pertcolvarres : conf.special().colvarres;
  if (!state.average_state_present) return true;
  size_t expected = 0;
  for (size_t i = 0; i < specs.size(); i += perturbed ? 2 : 1)
    if (specs[i].averaging != 0 && specs[i].tau > 0.0) ++expected;
  if (expected == state.average_states.size()) return true;
  io::messages.add("CV averaging restart count does not match the restraint definitions. Remove the CV state block only if resetting averages is intended.",
                   "Colvar", io::message::error);
  return false;
}

bool interaction::restore_colvar_average(Colvar_Bias &bias,
    const simulation::Parameter::colvar_bias_spec &spec, unsigned int ordinal,
    double timestep, const configuration::Configuration &conf, bool perturbed)
{
  const auto &state = perturbed ? conf.special().pertcolvarres : conf.special().colvarres;
  if (!state.average_state_present || !bias.settings().use_time_average) return true;
  for (const auto &record : state.average_states) {
    const bool identity = spec.index != 0 ? record.index == spec.index
                                         : record.index == 0 && record.ordinal == ordinal;
    if (!identity) continue;
    if (record.type != spec.type || record.method != spec.method
        || record.tau != spec.tau || record.timestep != timestep
        || !bias.restore_average_state(record.average, record.initialised)) {
      io::messages.add("CV averaging restart metadata/state mismatch for restraint row "
                       + std::to_string(ordinal) + ". Keep IDs, types, averaging settings and timestep consistent.",
                       "Colvar", io::message::error);
      return false;
    }
    return true;
  }
  io::messages.add("Missing CV averaging restart state for restraint row " + std::to_string(ordinal),
                   "Colvar", io::message::error);
  return false;
}

void interaction::checkpoint_colvar_average(const Colvar_Bias &bias,
    const simulation::Parameter::colvar_bias_spec &spec, unsigned int ordinal,
    double timestep, configuration::Configuration &conf, bool perturbed)
{
  if (!bias.settings().use_time_average) return;
  auto &state = perturbed ? conf.special().pertcolvarres : conf.special().colvarres;
  configuration::Colvar_Average_State record;
  record.ordinal = ordinal;
  record.index = spec.index;
  record.type = spec.type;
  record.method = spec.method;
  record.tau = spec.tau;
  record.timestep = timestep;
  record.average = bias.average_state();
  record.initialised = bias.average_initialised();
  state.average_states.push_back(record);
  state.average_state_present = true;
}
