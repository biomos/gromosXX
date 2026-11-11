/**
 * @file perturbed_colvar_restraint_interaction.cc
 * perturbed collective variable restraining
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"


#include "../../interaction/interaction_types.h"

#include "../../interaction/special/perturbed_colvar_restraint_interaction.h"
#include "../../interaction/special/colvar/angle.h"
#include "../../interaction/special/colvar/colvar_state.h"
#include "../../interaction/special/colvar/coordnum.h"
#include "../../interaction/special/colvar/dihedral.h"
#include "../../interaction/special/colvar/distance.h"
#include "../../interaction/special/colvar/distance_difference.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

namespace {

int rah_mode_from_rah(int rah)
{
  if (rah >= -1 && rah <= 1) return rah;
  if (rah >= 9 && rah <= 11) return rah - 10;
  if (rah >= 19 && rah <= 21) return rah - 20;
  if (rah >= 29 && rah <= 31) return rah - 30;
  if (rah >= 39 && rah <= 41) return rah - 40;
  if (rah >= 49 && rah <= 51) return rah - 50;
  if (rah >= 59 && rah <= 61) return rah - 60;
  return 0;
}

interaction::Colvar_Bias::Wall_Mode wall_mode_from_rah(int rah)
{
  const int rah_mode = rah_mode_from_rah(rah);
  if (rah_mode < 0) return interaction::Colvar_Bias::wall_lower;
  if (rah_mode > 0) return interaction::Colvar_Bias::wall_upper;
  return interaction::Colvar_Bias::wall_two_sided;
}

bool angular_type(const std::string &type)
{
  return type == "ANGLE" || type == "DIHEDRAL";
}

const simulation::Parameter::colvar_bias_spec *nth_spec(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type,
    size_t index)
{
  size_t n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->type != type) continue;
    if (n == index) return &(*it);
    ++n;
  }
  return NULL;
}

size_t count_specs(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type)
{
  size_t n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->type == type) ++n;
  }
  return n;
}

double interpolate(double A, double B, double lambda)
{
  return (1.0 - lambda) * A + lambda * B;
}


} // anonymous namespace

interaction::Perturbed_Colvar_Restraint_Interaction::~Perturbed_Colvar_Restraint_Interaction()
{
  for (std::vector<Term>::iterator it = m_terms.begin(), to = m_terms.end();
       it != to; ++it) {
    delete it->cv;
    delete it->bias;
  }
}

interaction::Colvar_Bias::Settings
interaction::Perturbed_Colvar_Restraint_Interaction::settings_from_spec(
  const simulation::Parameter::colvar_bias_spec &spec,
  simulation::Simulation &sim) const
{
  Colvar_Bias::Settings settings;

  settings.target = spec.target;
  settings.k = spec.k;
  settings.weight = 1.0;
  settings.linear_cutoff = spec.linear_tail;
  settings.rah = spec.rah;
  settings.virial = spec.virial;
  settings.wall_mode = wall_mode_from_rah(spec.rah);
  settings.periodic = (spec.type == "DIHEDRAL");
  settings.period = settings.periodic ? 2.0 * math::Pi : 0.0;
  settings.paired_wall = spec.paired_wall;
  settings.wall_partner_target = spec.wall_partner_target;
  settings.explicit_periodic_branch = spec.explicit_periodic_branch;
  settings.branch_delta = spec.branch_delta;

  settings.use_time_average = (spec.averaging != 0);
  if (spec.averaging == 2) {
    settings.average_transform = Colvar_Bias::average_inverse_cubic;
  }
  else if (spec.averaging == 1) {
    settings.average_transform = Colvar_Bias::average_identity;
  }
  else {
    settings.average_transform = Colvar_Bias::average_none;
  }

  if (settings.use_time_average && spec.tau > 0.0) {
    settings.exponential_term = std::exp(-sim.time_step_size() / spec.tau);
    if (spec.force_scale == 1) {
      settings.force_scale_mode = Colvar_Bias::forcescale_relaxation;
    }
    else if (spec.force_scale == 2) {
      settings.force_scale_mode = Colvar_Bias::forcescale_chain_rule;
    }
    else {
      settings.force_scale_mode = Colvar_Bias::forcescale_none;
    }
  }
  else {
    settings.use_time_average = false;
    settings.average_transform = Colvar_Bias::average_none;
    settings.exponential_term = 0.0;
    settings.force_scale_mode = Colvar_Bias::forcescale_none;
  }

  return settings;
}

int interaction::Perturbed_Colvar_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
                         configuration::Configuration &conf,
                         simulation::Simulation &sim)
{
  m_timer.start(sim);
  const bool advance_average = m_average_step != static_cast<long long>(sim.steps());
  conf.special().pertcolvarres.average_states.clear();
  conf.special().pertcolvarres.average_state_present = false;
  double Ctot = 0.0;
  double Etot = 0.0;

  conf.special().pertcolvarres.energies.clear();
  conf.special().pertcolvarres.values.clear();
  conf.special().pertcolvarres.types.clear();

  for (size_t i = 0; i < m_terms.size(); ++i) {
    Term &term = m_terms[i];
    Colvar *cv = term.cv;

    util::Virtual_Atom *firstatom = cv->atoms[0];
    const unsigned int eg = topo.atom_energy_group()[firstatom->atom(0)];
    const double lambda =
      topo.individual_lambda(simulation::colvarres_lambda)[eg][eg];
    const double lambda_derivative =
      topo.individual_lambda_derivative(simulation::colvarres_lambda)[eg][eg];

    if (cv->calculate(topo, conf, sim)) {
      std::cerr << "Colvar: error during calculation of "
                << cv->name << "\n";
      io::messages.add("Error in perturbed colvar calc",
                       "Colvar", io::message::error);
      m_timer.stop();
      return 1;
    }

    Ctot += cv->ct;
    conf.special().pertcolvarres.values.push_back(cv->ct);
    conf.special().pertcolvarres.types.push_back(term.type);

    if (!sum) {
      const double E = apply_restraint(topo, conf, sim, term,
                                       lambda, lambda_derivative, eg, advance_average);
      if (!std::isfinite(E)) {
        m_timer.stop();
        return 1;
      }
      conf.special().pertcolvarres.energies.push_back(E);
      Etot += E;
      conf.current().energies.add_colvar_energy(term.type, true, eg, E);
    }
  }

  if (sum) {
    io::messages.add("Perturbed COLVARRES sum restraints are not implemented.",
                     "Colvar", io::message::error);
  }

  for (const auto &term : m_terms)
    checkpoint_colvar_average(*term.bias, term.A, term.ordinal, sim.time_step_size(), conf, true);
  m_average_step = sim.steps();
  conf.special().pertcolvarres.totv = Ctot;
  conf.special().pertcolvarres.tote = Etot;
  m_timer.stop();
  return 0;
}

bool interaction::Perturbed_Colvar_Restraint_Interaction::evaluate_bias(
    const Term &term, Colvar_Bias &bias, simulation::Simulation &sim,
    double lambda, double lambda_derivative, bool advance_average,
    Bias_Evaluation &evaluation) const
{
  const auto envelope = Colvar_Bias::hidden_envelope(
    lambda, term.A.hidden_n, term.A.hidden_m);
  if (!std::isfinite(lambda) || !std::isfinite(lambda_derivative)
      || !std::isfinite(envelope.value) || !std::isfinite(envelope.derivative)) {
    io::messages.add("PERTCOLVARRES: non-finite lambda or hidden-restraint envelope at row "
                     + std::to_string(term.ordinal) + ". Check lambda mapping and N M.",
                     "Colvar", io::message::error);
    return false;
  }
  Colvar_Bias::Settings settings = settings_from_spec(term.A, sim);
  settings.target = interpolate(term.A.target, term.B.target, lambda);
  settings.wall_partner_target = interpolate(term.A.wall_partner_target,
                                             term.B.wall_partner_target, lambda);
  settings.k = interpolate(term.A.k, term.B.k, lambda);
  settings.linear_cutoff =
    interpolate(term.A.linear_tail, term.B.linear_tail, lambda);
  settings.virial = term.A.virial || term.B.virial;

  bias.update_settings_preserve_state(settings);
  Colvar_Bias::Result result = bias.evaluate(term.cv->ct, advance_average);

  const double D_target = term.B.target - term.A.target;
  const double D_k = term.B.k - term.A.k;
  const double dpotdl = result.dE_dk * D_k - result.dE_dvalue * D_target
    + result.dE_dlinear_cutoff * (term.B.linear_tail - term.A.linear_tail);
  const double energy_derivative = lambda_derivative
    * (envelope.derivative * result.energy + envelope.value * dpotdl);

  const double energy = envelope.value * result.energy;
  const double force_derivative = envelope.value * result.dE_dinstant;
  if (!std::isfinite(energy) || !std::isfinite(force_derivative)
      || !std::isfinite(energy_derivative)) {
    io::messages.add("PERTCOLVARRES: non-finite scaled bias at row "
                     + std::to_string(term.ordinal) + ". Check bias parameters and N M.",
                     "Colvar", io::message::error);
    return false;
  }
  evaluation.energy = energy;
  evaluation.derivative = energy_derivative;
  evaluation.force_derivative = force_derivative;
  return true;
}

bool interaction::Perturbed_Colvar_Restraint_Interaction::sample_lambdas(
    const topology::Topology &topo, configuration::Configuration &conf,
    simulation::Simulation &sim, const Term &term, unsigned int group) const
{
  const unsigned int count = sim.param().precalclam.nr_lambdas;
  if (count == 0) return true;
  if (count < 2 || conf.current().energies.AB_colvar.size() != count
      || conf.current().perturbed_energy_derivatives.AB_colvar.size() != count) {
    io::messages.add("COLVAR extended TI: invalid lambda grid or energy buffer size.",
                     "Colvar", io::message::error);
    return false;
  }
  for (unsigned int i = 0; i < count; ++i) {
    const double lambda = sim.param().precalclam.min_lam
      + (sim.param().precalclam.max_lam - sim.param().precalclam.min_lam)
        * double(i) / double(count - 1);
    const auto mapped = topo.individual_lambda_at(simulation::colvarres_lambda, group, group, lambda);
    // Reuse the live history after its one physical update. Counterfactual
    // samples must not change that history, live settings, topology or forces.
    Colvar_Bias sample(*term.bias);
    Bias_Evaluation evaluation;
    if (!evaluate_bias(term, sample, sim, mapped.first, mapped.second, false, evaluation)) return false;
    conf.current().energies.AB_colvar[i] += evaluation.energy;
    conf.current().perturbed_energy_derivatives.AB_colvar[i] += evaluation.derivative;
  }
  return true;
}

double interaction::Perturbed_Colvar_Restraint_Interaction::apply_restraint(
    topology::Topology &topo, configuration::Configuration &conf,
    simulation::Simulation &sim, Term &term, double lambda,
    double lambda_derivative, unsigned int energy_group, bool advance_average)
{
  Bias_Evaluation evaluation;
  if (!evaluate_bias(term, *term.bias, sim, lambda, lambda_derivative, advance_average, evaluation))
    return std::numeric_limits<double>::quiet_NaN();
  if (!sample_lambdas(topo, conf, sim, term, energy_group))
    return std::numeric_limits<double>::quiet_NaN();
  term.cv->apply_forces(topo, conf, evaluation.force_derivative,
    term.A.virial || term.B.virial, sim.param().pcouple.virial==math::molecular_virial);
  conf.current().perturbed_energy_derivatives.add_colvar_energy(
    term.type, true, energy_group, evaluation.derivative);
  return evaluation.energy;
}

int interaction::Perturbed_Colvar_Restraint_Interaction::init(
  topology::Topology &topo,
  configuration::Configuration &conf,
  simulation::Simulation &sim,
  std::ostream &os,
  bool quiet)
{
  m_terms.clear();
  m_average_step = -1;

  const std::vector<simulation::Parameter::colvar_bias_spec> &specs =
    sim.param().colvarres.pert_bias_specs;

  if (!check_colvar_average_count(conf, specs, true)) return 1;

  if (specs.empty()) {
    if (!quiet) {
      os << "Perturbed colvar restraint interaction (no PERTRESTRAINTS)"
         << std::endl;
    }
    return 0;
  }

  size_t distance_index =
    count_specs(sim.param().colvarres.bias_specs, "DISTANCE");
  size_t difference_index =
    count_specs(sim.param().colvarres.bias_specs, "DISTANCE_DIFFERENCE");
  size_t angle_index =
    count_specs(sim.param().colvarres.bias_specs, "ANGLE");
  size_t dihedral_index =
    count_specs(sim.param().colvarres.bias_specs, "DIHEDRAL");
  size_t coordnum_index =
    count_specs(sim.param().colvarres.bias_specs, "COORDNUM");

  const char *known_types[] = {"DISTANCE", "ANGLE", "DIHEDRAL", "COORDNUM", "DISTANCE_DIFFERENCE"};
  for (size_t t = 0; t < sizeof(known_types) / sizeof(known_types[0]); ++t) {
    const size_t n = count_specs(specs, known_types[t]);
    if (n % 2 != 0) {
      io::messages.add(std::string("PERTCOLVARRES: odd number of ")
                       + known_types[t]
                       + " RESTRAINTS entries; last entry will use identical A/B states.",
                       "Perturbed_Colvar_Restraint_Interaction",
                       io::message::warning);
    }
  }

  for (size_t i = 0; i < specs.size(); ++i) {
    const simulation::Parameter::colvar_bias_spec &A = specs[i];
    const simulation::Parameter::colvar_bias_spec *B = nth_spec(specs, A.type, 1);
    size_t type_instance = 0;
    for (size_t j = 0; j < i; ++j) {
      if (specs[j].type == A.type) ++type_instance;
    }
    if (type_instance % 2 != 0) continue;
    B = nth_spec(specs, A.type, type_instance + 1);

    Term term;
    term.A = A;
    term.ordinal = i / 2 + 1;
    if (B == NULL) {
      io::messages.add("PERTCOLVARRES: missing B-state RESTRAINTS entry for type "
                       + A.type + "; using A state as B state.",
                       "Perturbed_Colvar_Restraint_Interaction",
                       io::message::warning);
      term.B = A;
    }
    else {
      term.B = *B;
    }
    term.type = A.type;

    if (A.type == "DISTANCE") {
      if (distance_index >= topo.distance_restraints().size()) continue;
      Distance_Colvar *cv = new Distance_Colvar();
      cv->params = &topo.distance_restraints()[distance_index++];
      term.cv = cv;
    }
    else if (A.type == "DISTANCE_DIFFERENCE") {
      if (difference_index >= topo.distance_difference_restraints().size()) {
        io::messages.add("Missing PERTDISTANCEDIFFRESSPEC geometry.", "Colvar", io::message::error);
        return 1;
      }
      Distance_Difference_Colvar *cv = new Distance_Difference_Colvar();
      cv->params = &topo.distance_difference_restraints()[difference_index++];
      term.cv = cv;
    }
    else if (A.type == "ANGLE") {
      if (angle_index >= topo.colvar_angle_geometries().size()) continue;
      Angle_Colvar *cv = new Angle_Colvar();
      cv->params = &topo.colvar_angle_geometries()[angle_index++];
      term.cv = cv;
    }
    else if (A.type == "DIHEDRAL") {
      if (dihedral_index >= topo.colvar_dihedral_geometries().size()) continue;
      Dihedral_Colvar *cv = new Dihedral_Colvar();
      cv->params = &topo.colvar_dihedral_geometries()[dihedral_index++];
      term.cv = cv;
    }
    else if (A.type == "COORDNUM") {
      if (coordnum_index >= topo.coordnum_restraint().size()) continue;
      Coordnum_Colvar *cv = new Coordnum_Colvar();
      cv->params = &topo.coordnum_restraint()[coordnum_index++];
      term.cv = cv;
    }
    else {
      continue;
    }

    if (term.cv->init(topo, conf, sim, os, quiet)) {
      io::messages.add("Error in perturbed colvar init",
                       "Colvar", io::message::error);
      delete term.cv;
      return 1;
    }

    term.bias = new Colvar_Bias(settings_from_spec(term.A, sim));
    if (!restore_colvar_average(*term.bias, term.A, term.ordinal, sim.time_step_size(), conf, true)) {
      delete term.bias;
      delete term.cv;
      return 1;
    }

    DEBUG(9, "PERTCOLVARRES init cv " << m_terms.size() + 1
          << " index " << term.A.index
          << " type " << term.type
          << " method " << term.A.method
          << " targetA " << term.A.target
          << " targetB " << term.B.target
          << " KA " << term.A.k
          << " KB " << term.B.k);

    m_terms.push_back(term);
  }

  conf.special().pertcolvarres.average_states.clear();
  conf.special().pertcolvarres.average_state_present = false;
  for (const auto &term : m_terms)
    checkpoint_colvar_average(*term.bias, term.A, term.ordinal, sim.time_step_size(), conf, true);

  if (!quiet) {
    os << "Perturbed colvar restraint interaction" << std::endl;
  }
  return 0;
}
