/** EDS distance biases using the common CV geometry, bias and force layers. */
#include "../../stdheader.h"
#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../configuration/configuration.h"
#include "../../simulation/simulation.h"
#include "colvar/colvar.h"
#include "colvar/distance.h"
#include "colvar/colvar_bias.h"
#include "eds_colvar_restraint_interaction.h"

struct interaction::Eds_Colvar_Restraint_Interaction::Term {
  topology::distance_restraint_struct geometry;
  Distance_Colvar cv;
  std::vector<Colvar_Bias> biases;
  unsigned int index;
  bool virial;
  Term(const topology::eds_distance_restraint_struct &source)
    : geometry(source.v1,source.v2,source.r0[0],1.0,source.rah), index(0), virial(false) {
    cv.params = &geometry;
  }
};

interaction::Eds_Colvar_Restraint_Interaction::~Eds_Colvar_Restraint_Interaction()
{
  for (auto term : m_terms) delete term;
}

int interaction::Eds_Colvar_Restraint_Interaction::init(topology::Topology &topo,
    configuration::Configuration &conf, simulation::Simulation &sim,
    std::ostream &os, bool quiet)
{
  for (auto term : m_terms) delete term;
  m_terms.clear();
  auto &diagnostics = conf.special().colvarres;
  diagnostics.eds_indices.clear();
  diagnostics.eds_values.clear();
  diagnostics.eds_energies.clear();
  const auto &specs = sim.param().colvarres.eds_bias_specs;
  if (specs.empty()) return 0;
  const unsigned int states = sim.param().eds.numstates;
  if (!sim.param().eds.eds || sim.param().eds.form < simulation::single_s
      || sim.param().eds.form > simulation::aeds_advanced_search2
      || states < 2 || specs.size() != topo.eds_distance_restraints().size()) {
    io::messages.add("EDS CV definitions require enabled EDS and matching geometry/state counts.","Colvar",io::message::error);
    return 1;
  }
  for (size_t i = 0; i < specs.size(); ++i) {
    if (specs[i].size() != states || topo.eds_distance_restraints()[i].r0.size() != states) return 1;
    Term *term = new Term(topo.eds_distance_restraints()[i]);
    if (term->cv.init(topo,conf,sim,os,true)) { delete term; return 1; }
    term->index = specs[i][0].index;
    term->virial = specs[i][0].virial;
    for (const auto &spec : specs[i]) {
      if (spec.method != 0 || spec.tau != 0 || spec.type != "DISTANCE") {
        delete term;
        io::messages.add("EDS CV currently supports instantaneous DISTANCE only.","Colvar",io::message::error);
        return 1;
      }
      Colvar_Bias::Settings settings;
      settings.target = spec.target;
      settings.k = spec.k;
      settings.linear_cutoff = spec.linear_tail;
      settings.wall_mode = static_cast<Colvar_Bias::Wall_Mode>(term->cv.rah_mode_value());
      term->biases.push_back(Colvar_Bias(settings));
    }
    m_terms.push_back(term);
  }
  if (!quiet) os << "EDS CV distance restraints: " << m_terms.size() << " terms, " << states << " states\n";
  return 0;
}

int interaction::Eds_Colvar_Restraint_Interaction::calculate_interactions(
    topology::Topology &topo, configuration::Configuration &conf, simulation::Simulation &sim)
{
  m_timer.start(sim);
  auto &diagnostics = conf.special().colvarres;
  diagnostics.eds_indices.clear();
  diagnostics.eds_values.clear();
  diagnostics.eds_energies.clear();
  const size_t states = sim.param().eds.numstates;
  if (!m_terms.empty() && (conf.current().energies.eds_vi.size() != states
      || conf.current().energies.eds_vi_special.size() != states
      || conf.special().eds.force_endstates.size() != states
      || conf.special().eds.virial_tensor_endstates.size() != states)) {
    io::messages.add("EDS CV state buffers have inconsistent sizes.","Colvar",io::message::error);
    m_timer.stop(); return 1;
  }
  if (!m_terms.empty()) {
    for (const auto &forces : conf.special().eds.force_endstates)
      if (forces.size() != topo.num_atoms()) {
        io::messages.add("EDS CV force buffer does not match the atom count.","Colvar",io::message::error);
        m_timer.stop(); return 1;
      }
  }
  for (auto term : m_terms) {
    if (term->cv.calculate(topo,conf,sim)) { m_timer.stop(); return 1; }
    std::vector<double> energies;
    for (size_t state = 0; state < states; ++state) {
      const auto bias = term->biases[state].evaluate(term->cv.ct);
      if (!std::isfinite(bias.energy) || !std::isfinite(bias.dE_dinstant)) {
        io::messages.add("Non-finite EDS CV bias.","Colvar",io::message::error);
        m_timer.stop(); return 1;
      }
      term->cv.apply_forces(topo,conf,bias.dE_dinstant,term->virial,
          conf.special().eds.force_endstates[state],conf.special().eds.virial_tensor_endstates[state],
          sim.param().pcouple.virial==math::molecular_virial);
      conf.current().energies.eds_vi[state] += bias.energy;
      conf.current().energies.eds_vi_special[state] += bias.energy;
      energies.push_back(bias.energy);
    }
    diagnostics.eds_indices.push_back(term->index);
    diagnostics.eds_values.push_back(term->cv.ct);
    diagnostics.eds_energies.push_back(energies);
  }
  m_timer.stop();
  return 0;
}
