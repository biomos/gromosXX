/**
 * @file colvar_restraint_interaction.cc
 * template methods of Colvar_Restraint_Interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"


// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/colvar_restraint_interaction.h"
#include "../../interaction/special/colvar/colvar.h"
#include "../../interaction/special/colvar/colvar_bias.h"
#include "../../interaction/special/colvar/colvar_state.h"
#include "../../interaction/special/colvar/angle.h"
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
  // Legacy DISTANCERESSPEC stores full/repulsive/attractive behavior in RAH.
  // Preserve that convention when distance geometries are reused as colvars.
  const int rah_mode = rah_mode_from_rah(rah);
  if (rah_mode < 0) {
    return interaction::Colvar_Bias::wall_lower;
  }
  if (rah_mode > 0) {
    return interaction::Colvar_Bias::wall_upper;
  }
  return interaction::Colvar_Bias::wall_two_sided;
}

const char *wall_mode_name(interaction::Colvar_Bias::Wall_Mode mode)
{
  if (mode == interaction::Colvar_Bias::wall_lower) {
    return "half-harmonic-repulsive";
  }
  if (mode == interaction::Colvar_Bias::wall_upper) {
    return "half-harmonic-attractive";
  }
  return "full-harmonic";
}

const char *average_transform_name(interaction::Colvar_Bias::Average_Transform transform)
{
  if (transform == interaction::Colvar_Bias::average_identity) {
    return "time";
  }
  if (transform == interaction::Colvar_Bias::average_inverse_cubic) {
    return "inverse-cubic";
  }
  return "none";
}

const char *force_scale_name(interaction::Colvar_Bias::Force_Scale_Mode mode)
{
  if (mode == interaction::Colvar_Bias::forcescale_relaxation) {
    return "relaxation";
  }
  if (mode == interaction::Colvar_Bias::forcescale_chain_rule) {
    return "chain-rule";
  }
  return "none";
}

bool angular_type(const std::string &type)
{
  return type == "ANGLE" || type == "DIHEDRAL";
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


} // anonymous namespace

interaction::Colvar_Restraint_Interaction::~Colvar_Restraint_Interaction(){
  for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
        delete *it;
  }
  for (std::vector<Colvar_Bias *>::iterator it = m_biases.begin(),
        to = m_biases.end(); it != to; ++it) {
        delete *it;
  }
}

interaction::Colvar_Bias::Settings
interaction::Colvar_Restraint_Interaction::settings_from_spec(
  const simulation::Parameter::colvar_bias_spec &spec,
  simulation::Simulation &sim) const
{
  Colvar_Bias::Settings settings;

  settings.target = spec.target;
  settings.k = spec.k;
  settings.weight = 1.0;
  settings.rah = spec.rah;
  settings.virial = spec.virial;

  // In the new COLVARRES format, LINEAR_TAIL is interpreted as the
  // harmonic displacement where the potential switches to a linear tail.
  // A value <= 0 disables linearisation.
  settings.linear_cutoff = spec.linear_tail;

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

  if (settings.use_time_average) {
    if (spec.tau <= 0.0) {
      io::messages.add("COLVARRES: averaging requested with TAU <= 0.0; disabling averaging for this CV.",
                       "Colvar", io::message::warning);
      settings.use_time_average = false;
      settings.average_transform = Colvar_Bias::average_none;
      settings.exponential_term = 0.0;
      settings.force_scale_mode = Colvar_Bias::forcescale_none;
    }
    else {
      settings.exponential_term = std::exp(- sim.time_step_size() / spec.tau);

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
  }
  else {
    settings.exponential_term = 0.0;
    settings.force_scale_mode = Colvar_Bias::forcescale_none;
  }

  return settings;
}

/**
 * calculate colvar restraint interactions
 */
int interaction::Colvar_Restraint_Interaction
::calculate_interactions(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  m_timer.start(sim);
  const bool advance_average = m_average_step != static_cast<long long>(sim.steps());
  conf.special().colvarres.average_states.clear();
  conf.special().colvarres.average_state_present = false;
  double Ctot=0;
  double Etot=0;
  conf.special().colvarres.energies.clear();
  conf.special().colvarres.values.clear();
  conf.special().colvarres.types.clear();

  // loop over the colvars and calculate values and derivatives
  for (size_t i = 0; i < m_colvars.size(); ++i) {
    Colvar *cv = m_colvars[i];
    Colvar_Bias *bias = m_biases[i];

    m_timer.start_subtimer("calculate");
    if (cv->calculate(topo, conf, sim)) {
       std::cerr << "Colvar: error during calculation of "
         << cv->name << "\n";

      io::messages.add("Error in colvar calc",
                       "Colvar",
                       io::message::error);
      m_timer.stop_subtimer("calculate");
      m_timer.stop();
      return 1;
    }
    m_timer.stop_subtimer("calculate");

    Ctot += cv->ct;
    conf.special().colvarres.values.push_back(cv->ct);
    if (i < m_colvar_types.size()) {
      conf.special().colvarres.types.push_back(m_colvar_types[i]);
    }
    else {
      conf.special().colvarres.types.push_back(cv->name);
    }

    if (!sum) {
      double E = apply_restraint(topo, conf, sim,
                                 *cv,
                                 i < m_colvar_types.size() ? m_colvar_types[i] : cv->name,
                                 *bias, advance_average);
      conf.special().colvarres.energies.push_back(E);
      Etot += E;

      // add to the energy group of the first atom of first list
      util::Virtual_Atom *firstatom = (cv->atoms)[0];
      conf.current().energies.add_colvar_energy(m_colvar_types[i], false,
        topo.atom_energy_group()[firstatom->atom(0)], E);
    }
  }

  conf.special().colvarres.totv=Ctot;

  // Restraining the sum is intentionally not wired to the new per-CV
  // COLVARRES format yet, because there is no single unambiguous per-term K.
  if (sum) {
    io::messages.add("COLVARRES sum restraints are not implemented for the new per-CV bias format.",
                     "Colvar", io::message::error);
  }

  const auto &specs = sim.param().colvarres.bias_specs;
  for (size_t i = 0; i < m_biases.size(); ++i)
    checkpoint_colvar_average(*m_biases[i], specs[m_spec_ordinals[i]-1], m_spec_ordinals[i],
                              sim.time_step_size(), conf, false);
  m_average_step = sim.steps();
  conf.special().colvarres.tote=Etot;
  m_timer.stop();
  return 0;
}


// apply forces and calculate potential using the bias object assigned to this colvar
double interaction::Colvar_Restraint_Interaction
  ::apply_restraint(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 const Colvar &cv,
 const std::string &type,
 Colvar_Bias &bias, bool advance_average) {

  const double curr = cv.ct;

  double E = 0.0;

  if (sim.param().colvarres.colvarres == simulation::colvar_restr_off) {

       std::cerr << "Colvar restraints are off, yet you ended up in colvar's apply_restraint .. something is wrong ..\n";

      io::messages.add("Error in colvar calc",
                       "Colvar",
                       io::message::error);
      return E;
  }

  Colvar_Bias::Result result = bias.evaluate(curr, advance_average);
  E = result.energy;


  const Colvar_Bias::Settings &settings = bias.settings();
  DEBUG(9, "COLVARRES step " << sim.steps()
        << " value " << curr
        << " restrained_value " << result.value
        << " target " << settings.target
          << " rah " << settings.rah
          << " form " << wall_mode_name(settings.wall_mode)
          << " active " << result.active
        << " energy " << result.energy
        << " linear " << result.linear
        << " dE/dinstant " << result.dE_dinstant
        << " force_scale " << result.force_scale);

  if (angular_type(type)) {
    DEBUG(9, "COLVARRES angular step " << sim.steps()
          << " type " << type
          << " value_deg " << curr * 180.0 / math::Pi
          << " restrained_value_deg " << result.value * 180.0 / math::Pi
          << " target_deg " << settings.target * 180.0 / math::Pi);
  }

  if (result.linear) {
    DEBUG(9, "COLVARRES linearized step " << sim.steps()
          << " type " << type
          << " value " << curr
          << " restrained_value " << result.value
          << " target " << settings.target
          << " abs_deviation " << std::fabs(result.value - settings.target)
          << " cvlin " << settings.linear_cutoff
          << " energy " << result.energy
          << " dE/dinstant " << result.dE_dinstant);
  }

  cv.apply_forces(topo, conf, result.dE_dinstant, settings.virial,
                  sim.param().pcouple.virial==math::molecular_virial);

  return E;
}

int interaction::Colvar_Restraint_Interaction::init(topology::Topology &topo,
             configuration::Configuration &conf,
             simulation::Simulation &sim,
             std::ostream &os,
             bool quiet)
{
  m_average_step = -1;
  m_spec_ordinals.clear();
  size_t distance_index = 0;
  size_t difference_index = 0;
  size_t angle_index = 0;
  size_t dihedral_index = 0;
  size_t coordnum_index = 0;
  m_colvar_types.clear();

  const std::vector<simulation::Parameter::colvar_bias_spec> &specs =
    sim.param().colvarres.bias_specs;
  if (!check_colvar_average_count(conf, specs, false)) return 1;
  const size_t pert_distance_count =
    count_specs(sim.param().colvarres.pert_bias_specs, "DISTANCE") / 2;
  const size_t pert_angle_count =
    count_specs(sim.param().colvarres.pert_bias_specs, "ANGLE") / 2;
  const size_t pert_dihedral_count =
    count_specs(sim.param().colvarres.pert_bias_specs, "DIHEDRAL") / 2;
  const size_t pert_coordnum_count =
    count_specs(sim.param().colvarres.pert_bias_specs, "COORDNUM") / 2;

  for (size_t i = 0; i < specs.size(); ++i) {
    const simulation::Parameter::colvar_bias_spec &spec = specs[i];

    Colvar *cv = NULL;

    if (spec.type == "DISTANCE") {
      if (distance_index >= topo.distance_restraints().size()) {
        io::messages.add("COLVARRES contains more DISTANCE entries than DISTANCERESSPEC geometries.",
                         "Colvar", io::message::error);
        continue;
      }

      interaction::Distance_Colvar *distance_cv = new interaction::Distance_Colvar();
      distance_cv->params = &topo.distance_restraints()[distance_index++];
      cv = distance_cv;
    }
    else if (spec.type == "DISTANCE_DIFFERENCE") {
      if (difference_index >= topo.distance_difference_restraints().size()) {
        io::messages.add("Missing DISTANCEDIFFRESSPEC geometry.", "Colvar", io::message::error);
        return 1;
      }
      Distance_Difference_Colvar *difference_cv = new Distance_Difference_Colvar();
      difference_cv->params = &topo.distance_difference_restraints()[difference_index++];
      cv = difference_cv;
    }
    else if (spec.type == "ANGLE") {
      if (angle_index >= topo.colvar_angle_geometries().size()) {
        io::messages.add("COLVARRES contains more ANGLE entries than ANGRESSPEC geometries.",
                         "Colvar", io::message::error);
        continue;
      }

      interaction::Angle_Colvar *angle_cv = new interaction::Angle_Colvar();
      angle_cv->params = &topo.colvar_angle_geometries()[angle_index++];
      cv = angle_cv;
    }
    else if (spec.type == "DIHEDRAL") {
      if (dihedral_index >= topo.colvar_dihedral_geometries().size()) {
        io::messages.add("COLVARRES contains more DIHEDRAL entries than DIHEDRALRESSPEC geometries.",
                         "Colvar", io::message::error);
        continue;
      }

      interaction::Dihedral_Colvar *dihedral_cv = new interaction::Dihedral_Colvar();
      dihedral_cv->params = &topo.colvar_dihedral_geometries()[dihedral_index++];
      cv = dihedral_cv;
    }
    else if (spec.type == "COORDNUM") {
      if (coordnum_index >= topo.coordnum_restraint().size()) {
        io::messages.add("COLVARRES contains more COORDNUM entries than COORDNUMRESSPEC geometries.",
                         "Colvar", io::message::error);
        continue;
      }

      interaction::Coordnum_Colvar *coordnum_cv = new interaction::Coordnum_Colvar();
      coordnum_cv->params = &topo.coordnum_restraint()[coordnum_index++];
      cv = coordnum_cv;
    }
    else {
      io::messages.add("COLVARRES: unknown TYPE '" + spec.type + "'. Use DISTANCE, COORDNUM, ANGLE, DIHEDRAL, or DISTANCE_DIFFERENCE.",
                       "Colvar", io::message::error);
      continue;
    }

    if (cv->init(topo, conf, sim, os, quiet)) {
       os << "Colvar: error during initialisation of "
          << cv->name << "\n";

      io::messages.add("Error in colvar init",
                       "Colvar",
                       io::message::error);
      delete cv;
      return 1;
    }

    Ctot0 += spec.target;

    Colvar_Bias::Settings settings = settings_from_spec(spec, sim);

    DEBUG(9, "COLVARRES init cv " << i + 1
          << " index " << spec.index
          << " type " << spec.type
          << " method " << spec.method
          << " target " << settings.target
          << " K " << settings.k
          << " weight " << settings.weight
          << " rah " << spec.rah
          << " form " << wall_mode_name(settings.wall_mode)
          << " linear_cutoff " << settings.linear_cutoff
          << " virial " << settings.virial
          << " averaging " << average_transform_name(settings.average_transform)
          << " use_time_average " << settings.use_time_average
          << " tau_exp " << settings.exponential_term
          << " force_scale " << force_scale_name(settings.force_scale_mode));

    if (angular_type(spec.type)) {
      DEBUG(9, "COLVARRES init angular cv " << i + 1
            << " type " << spec.type
            << " target_deg " << settings.target * 180.0 / math::Pi);
    }

    Colvar_Bias *bias = new Colvar_Bias(settings);

    if (settings.use_time_average && conf.special().colvarres.average_state_present) {
      if (!restore_colvar_average(*bias, spec, i+1, sim.time_step_size(), conf, false)) {
        delete bias;
        delete cv;
        return 1;
      }
    } else if (settings.use_time_average) {
      // Initialise the average from the current instantaneous value.
      if (cv->calculate(topo, conf, sim)) {
        io::messages.add("Error initialising colvar average",
                         "Colvar", io::message::error);
        delete bias;
        delete cv;
        return 1;
      }
      bias->reset_average(cv->ct);
    }

    m_colvars.push_back(cv);
    m_biases.push_back(bias);
    m_spec_ordinals.push_back(i+1);
    m_colvar_types.push_back(spec.type);
  }

  conf.special().colvarres.average_states.clear();
  conf.special().colvarres.average_state_present = false;
  for (size_t i = 0; i < m_biases.size(); ++i)
    checkpoint_colvar_average(*m_biases[i], specs[m_spec_ordinals[i]-1], m_spec_ordinals[i],
                              sim.time_step_size(), conf, false);

  if (distance_index + pert_distance_count < topo.distance_restraints().size()) {
    io::messages.add("DISTANCERESSPEC contains more geometries than COLVARRES DISTANCE entries; extra geometries are ignored.",
                     "Colvar", io::message::warning);
  }

  if (angle_index + pert_angle_count < topo.colvar_angle_geometries().size()) {
    io::messages.add("ANGRESSPEC contains more geometries than COLVARRES ANGLE entries; extra geometries are ignored.",
                     "Colvar", io::message::warning);
  }

  if (dihedral_index + pert_dihedral_count < topo.colvar_dihedral_geometries().size()) {
    io::messages.add("DIHEDRALRESSPEC contains more geometries than COLVARRES DIHEDRAL entries; extra geometries are ignored.",
                     "Colvar", io::message::warning);
  }

  if (coordnum_index + pert_coordnum_count < topo.coordnum_restraint().size()) {
    io::messages.add("COORDNUMRESSPEC contains more geometries than COLVARRES COORDNUM entries; extra geometries are ignored.",
                     "Colvar", io::message::warning);
  }

  if (!quiet) {
    os << "Colvar restraint interaction";
    os << std::endl;
  }
  return 0;
}
