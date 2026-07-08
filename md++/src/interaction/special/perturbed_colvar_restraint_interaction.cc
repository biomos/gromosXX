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

#include "../../math/periodicity.h"

#include "../../interaction/interaction_types.h"

#include "../../interaction/special/perturbed_colvar_restraint_interaction.h"
#include "../../interaction/special/colvar/angle.h"
#include "../../interaction/special/colvar/coordnum.h"
#include "../../interaction/special/colvar/dihedral.h"
#include "../../interaction/special/colvar/distance.h"

#include "../../util/template_split.h"
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

template<math::boundary_enum B>
static void _apply_colvar_virial(
  topology::Topology &topo,
  configuration::Configuration &conf,
  const std::vector<util::Virtual_Atom*> &atoms,
  const std::vector<math::Vec> &forces)
{
  if (atoms.size() < 2 || forces.size() != atoms.size()) return;

  math::Periodicity<B> periodicity(conf.current().box);
  const math::Vec ref = atoms[0]->pos(conf, topo);

  for (size_t i = 1; i < atoms.size(); ++i) {
    math::Vec rel;
    periodicity.nearest_image(atoms[i]->pos(conf, topo), ref, rel);

    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        conf.current().virial_tensor(a, b) += rel(a) * forces[i](b);
      }
    }
  }
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
    }

    Ctot += cv->ct;
    conf.special().pertcolvarres.values.push_back(cv->ct);
    conf.special().pertcolvarres.types.push_back(term.type);

    if (!sum) {
      const double E = apply_restraint(topo, conf, sim, term,
                                       lambda, lambda_derivative, eg);
      conf.special().pertcolvarres.energies.push_back(E);
      Etot += E;
      conf.current().energies.colvarres_energy[eg] += E;
    }
  }

  if (sum) {
    io::messages.add("Perturbed COLVARRES sum restraints are not implemented.",
                     "Colvar", io::message::error);
  }

  conf.special().pertcolvarres.totv = Ctot;
  conf.special().pertcolvarres.tote = Etot;
  m_timer.stop();
  return 0;
}

double interaction::Perturbed_Colvar_Restraint_Interaction
::apply_restraint(topology::Topology &topo,
                  configuration::Configuration &conf,
                  simulation::Simulation &sim,
                  Term &term,
                  double lambda,
                  double lambda_derivative,
                  unsigned int energy_group)
{
  Colvar_Bias::Settings settings = settings_from_spec(term.A, sim);
  settings.target = interpolate(term.A.target, term.B.target, lambda);
  settings.k = interpolate(term.A.k, term.B.k, lambda);
  settings.linear_cutoff =
    interpolate(term.A.linear_tail, term.B.linear_tail, lambda);
  settings.virial = term.A.virial || term.B.virial;

  term.bias->update_settings_preserve_state(settings);
  Colvar_Bias::Result result = term.bias->evaluate(term.cv->ct);

  std::vector<math::Vec> forces(term.cv->atoms.size(), math::Vec(0));
  for (size_t i = 0; i < term.cv->atoms.size(); ++i) {
    forces[i] = -result.dE_dinstant * term.cv->derivatives[i];
    term.cv->atoms[i]->force(conf, topo, forces[i]);
  }

  if (settings.virial) {
    SPLIT_BOUNDARY(_apply_colvar_virial, topo, conf, term.cv->atoms, forces);
  }

  const double D_target = term.B.target - term.A.target;
  const double D_k = term.B.k - term.A.k;
  double dE_dk = 0.0;
  if (settings.k != 0.0) {
    dE_dk = result.energy / settings.k;
  }
  const double dpotdl = dE_dk * D_k - result.dE_dvalue * D_target;
  const double energy_derivative = lambda_derivative * dpotdl;

  conf.current().perturbed_energy_derivatives.
    colvarres_energy[energy_group] += energy_derivative;

  DEBUG(9, "PERTCOLVARRES step " << sim.steps()
        << " type " << term.type
        << " lambda " << lambda
        << " value " << term.cv->ct
        << " target " << settings.target
        << " K " << settings.k
        << " energy " << result.energy
        << " dE/dlambda " << energy_derivative);

  if (angular_type(term.type)) {
    DEBUG(9, "PERTCOLVARRES angular step " << sim.steps()
          << " type " << term.type
          << " value_deg " << term.cv->ct * 180.0 / math::Pi
          << " target_deg " << settings.target * 180.0 / math::Pi);
  }

  return result.energy;
}

int interaction::Perturbed_Colvar_Restraint_Interaction::init(
  topology::Topology &topo,
  configuration::Configuration &conf,
  simulation::Simulation &sim,
  std::ostream &os,
  bool quiet)
{
  m_terms.clear();

  const std::vector<simulation::Parameter::colvar_bias_spec> &specs =
    sim.param().colvarres.pert_bias_specs;

  if (specs.empty()) {
    if (!quiet) {
      os << "Perturbed colvar restraint interaction (no PERTRESTRAINTS)"
         << std::endl;
    }
    return 0;
  }

  size_t distance_index =
    count_specs(sim.param().colvarres.bias_specs, "DISTANCE");
  size_t angle_index =
    count_specs(sim.param().colvarres.bias_specs, "ANGLE");
  size_t dihedral_index =
    count_specs(sim.param().colvarres.bias_specs, "DIHEDRAL");
  size_t coordnum_index =
    count_specs(sim.param().colvarres.bias_specs, "COORDNUM");

  const char *known_types[] = {"DISTANCE", "ANGLE", "DIHEDRAL", "COORDNUM"};
  for (size_t t = 0; t < 4; ++t) {
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
    else if (A.type == "ANGLE") {
      if (angle_index >= topo.angle_restraints().size()) continue;
      Angle_Colvar *cv = new Angle_Colvar();
      cv->params = &topo.angle_restraints()[angle_index++];
      term.cv = cv;
    }
    else if (A.type == "DIHEDRAL") {
      if (dihedral_index >= topo.dihedral_restraints().size()) continue;
      Dihedral_Colvar *cv = new Dihedral_Colvar();
      cv->params = &topo.dihedral_restraints()[dihedral_index++];
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
    }

    term.bias = new Colvar_Bias(settings_from_spec(term.A, sim));

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

  if (!quiet) {
    os << "Perturbed colvar restraint interaction" << std::endl;
  }
  return 0;
}
