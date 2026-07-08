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

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/colvar_restraint_interaction.h"
#include "../../interaction/special/colvar/colvar.h"
#include "../../interaction/special/colvar/colvar_bias.h"
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

template<math::boundary_enum B>
static void _apply_colvar_virial(
  topology::Topology &topo,
  configuration::Configuration &conf,
  const std::vector<util::Virtual_Atom*> &atoms,
  const std::vector<math::Vec> &forces)
{
  if (atoms.size() < 2 || forces.size() != atoms.size()) {
    return;
  }

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
                                 cv->atoms,
                                 cv->derivatives,
                                 cv->ct,
                                 i < m_colvar_types.size() ? m_colvar_types[i] : cv->name,
                                 *bias);
      conf.special().colvarres.energies.push_back(E);
      Etot += E;

      // add to the energy group of the first atom of first list
      util::Virtual_Atom *firstatom = (cv->atoms)[0];
      conf.current().energies.colvarres_energy[topo.atom_energy_group()
        [(*firstatom).atom(0)]] += E;
    }
  }

  conf.special().colvarres.totv=Ctot;

  // Restraining the sum is intentionally not wired to the new per-CV
  // COLVARRES format yet, because there is no single unambiguous per-term K.
  if (sum) {
    io::messages.add("COLVARRES sum restraints are not implemented for the new per-CV bias format.",
                     "Colvar", io::message::error);
  }

  conf.special().colvarres.tote=Etot;
  m_timer.stop();
  return 0;
}


// apply forces and calculate potential using the bias object assigned to this colvar
double interaction::Colvar_Restraint_Interaction
  ::apply_restraint(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 const std::vector< util::Virtual_Atom* > &atoms,
 const math::VArray &derivatives,
 double curr,
 const std::string &type,
 Colvar_Bias &bias) {

  double E = 0.0;

  if (sim.param().colvarres.colvarres == simulation::colvar_restr_off) {

       std::cerr << "Colvar restraints are off, yet you ended up in colvar's apply_restraint .. something is wrong ..\n";

      io::messages.add("Error in colvar calc",
                       "Colvar",
                       io::message::error);
      return E;
  }

  Colvar_Bias::Result result = bias.evaluate(curr);
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

  std::vector<math::Vec> forces(atoms.size(), math::Vec(0));
  for (size_t i=0; i<atoms.size(); i++) {
    forces[i] = -result.dE_dinstant * derivatives[i];
    DEBUG(15, "COLVARRES step " << sim.steps()
          << " atom " << i + 1
          << " derivative " << math::v2s(derivatives[i])
          << " force " << math::v2s(forces[i]));
    (*atoms[i]).force(conf, topo, forces[i]);
  }

  if (settings.virial) {
    SPLIT_BOUNDARY(_apply_colvar_virial, topo, conf, atoms, forces);
  }

  return E;
}

int interaction::Colvar_Restraint_Interaction::init(topology::Topology &topo,
             configuration::Configuration &conf,
             simulation::Simulation &sim,
             std::ostream &os,
             bool quiet)
{
  size_t distance_index = 0;
  size_t angle_index = 0;
  size_t dihedral_index = 0;
  size_t coordnum_index = 0;
  m_colvar_types.clear();

  const std::vector<simulation::Parameter::colvar_bias_spec> &specs =
    sim.param().colvarres.bias_specs;
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
    else if (spec.type == "ANGLE") {
      if (angle_index >= topo.angle_restraints().size()) {
        io::messages.add("COLVARRES contains more ANGLE entries than ANGRESSPEC geometries.",
                         "Colvar", io::message::error);
        continue;
      }

      interaction::Angle_Colvar *angle_cv = new interaction::Angle_Colvar();
      angle_cv->params = &topo.angle_restraints()[angle_index++];
      cv = angle_cv;
    }
    else if (spec.type == "DIHEDRAL") {
      if (dihedral_index >= topo.dihedral_restraints().size()) {
        io::messages.add("COLVARRES contains more DIHEDRAL entries than DIHEDRALRESSPEC geometries.",
                         "Colvar", io::message::error);
        continue;
      }

      interaction::Dihedral_Colvar *dihedral_cv = new interaction::Dihedral_Colvar();
      dihedral_cv->params = &topo.dihedral_restraints()[dihedral_index++];
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
      io::messages.add("COLVARRES: unknown TYPE '" + spec.type + "'. Use DISTANCE, COORDNUM, ANGLE, or DIHEDRAL.",
                       "Colvar", io::message::error);
      continue;
    }

    if (cv->init(topo, conf, sim, os, quiet)) {
       os << "Colvar: error during initialisation of "
          << cv->name << "\n";

      io::messages.add("Error in colvar init",
                       "Colvar",
                       io::message::error);
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

    if (settings.use_time_average) {
      // Initialise the average from the current instantaneous value.
      cv->calculate(topo, conf, sim);
      bias->reset_average(cv->ct);
    }

    m_colvars.push_back(cv);
    m_biases.push_back(bias);
    m_colvar_types.push_back(spec.type);
  }

  if (distance_index + pert_distance_count < topo.distance_restraints().size()) {
    io::messages.add("DISTANCERESSPEC contains more geometries than COLVARRES DISTANCE entries; extra geometries are ignored.",
                     "Colvar", io::message::warning);
  }

  if (angle_index + pert_angle_count < topo.angle_restraints().size()) {
    io::messages.add("ANGRESSPEC contains more geometries than COLVARRES ANGLE entries; extra geometries are ignored.",
                     "Colvar", io::message::warning);
  }

  if (dihedral_index + pert_dihedral_count < topo.dihedral_restraints().size()) {
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
