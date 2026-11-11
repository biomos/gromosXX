/** Shared configuration checkpoint helpers for CV averaging histories. */
#ifndef INCLUDED_COLVAR_STATE_H
#define INCLUDED_COLVAR_STATE_H
#include "../../../simulation/multibath.h"
#include "../../../simulation/parameter.h"

namespace configuration { class Configuration; }
namespace interaction {
class Colvar_Bias;
bool check_colvar_average_count(const configuration::Configuration &conf,
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs, bool perturbed);
bool restore_colvar_average(Colvar_Bias &bias,
    const simulation::Parameter::colvar_bias_spec &spec, unsigned int ordinal,
    double timestep, const configuration::Configuration &conf, bool perturbed);
void checkpoint_colvar_average(const Colvar_Bias &bias,
    const simulation::Parameter::colvar_bias_spec &spec, unsigned int ordinal,
    double timestep, configuration::Configuration &conf, bool perturbed);
}
#endif
