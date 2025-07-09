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
#include "../../interaction/special/colvar/contactnum.h"
#include "../../interaction/special/colvar/distance.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

interaction::Colvar_Restraint_Interaction::~Colvar_Restraint_Interaction(){
  for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
        delete *it;
  }
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
  
  // loop over the colvars and calculate values and derivatives
  for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
    m_timer.start_subtimer("calculate");
    if ((*it)->calculate(topo, conf, sim)) {
       std::cerr << "Colvar: error during calculation of "
         << (*it)->name << "\n";

      io::messages.add("Error in colvar calc",
                       "Colvar",
                       io::message::error);
    }
    m_timer.stop_subtimer("calculate");
    Ctot+=(*it)->ct;
    conf.special().colvarres.values.push_back((*it)->ct);
    if (!sum) {
      double E;
      E= apply_restraint(topo, conf, sim, (*it)->atoms, (*it)->derivatives, (*it)->ct, (*it)->targetvalue, (*it)->w0, (*it)->rah);
      conf.special().colvarres.energies.push_back(E);
      Etot+=E;
      // add to the energy group of the first atom of first list
      util::Virtual_Atom *firstatom = ((*it)->atoms)[0];
      conf.current().energies.colvarres_energy[topo.atom_energy_group()
					    [(*firstatom).atom(0)]] += E;
    }
  }
  
  conf.special().colvarres.totv=Ctot;
  
  // not implemented yet, restrain to the sum of all colvar restraints
  if (sum) {
    for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
      double E;
      // the individual weights are ignored when restraining to the sum
      E= apply_restraint(topo, conf, sim, (*it)->atoms, (*it)->derivatives, Ctot, Ctot0, 1.0);
      conf.special().colvarres.energies.push_back(E);
      Etot+=E;
      // add to the energy group of the first atom of first list
      util::Virtual_Atom *firstatom = ((*it)->atoms)[0];
      conf.current().energies.colvarres_energy[topo.atom_energy_group()
					    [(*firstatom).atom(0)]] += E;
    }
  }
  conf.special().colvarres.tote=Etot;
  m_timer.stop();
  return 0;
}

// apply forces and calculate potential using the biasing function specified in the COLVARRES block
double interaction::Colvar_Restraint_Interaction
  ::apply_restraint(topology::Topology & topo,
 configuration::Configuration & conf,simulation::Simulation & sim, std::vector< util::Virtual_Atom* > atoms, 
                    math::VArray &derivatives, double &curr, double &target, double w0, int rah) {
  
  double K = sim.param().colvarres.K;
  double r_linear = sim.param().distanceres.r_linear;
  double diff = target - curr;
  double E = 0;

  DEBUG(8, "target value: " << target << ", current value: " << curr << ", rah: " << rah);
                  
  if (sim.param().colvarres.colvarres == simulation::colvar_restr_off) {   
       std::cerr << "Colvar restraints are off, yet you ended up in colvar's apply_restraint .. something is wrong ..\n";

      io::messages.add("Error in colvar calc",
                       "Colvar",
                       io::message::error);
  }
  // entering harmonic restraining
  else if (sim.param().colvarres.colvarres == simulation::colvar_restr_harmonic) {
  E=0.5*sim.param().colvarres.K*w0*diff*diff;
  DEBUG(8, "Harmonic restr energy: " << E);

  DEBUG(16, "Number of derivatives: " << derivatives.size());
  for (int i = 0; i < derivatives.size(); ++i) {
    DEBUG(8, "Derivative " << i << ": " << math::v2s(derivatives[i]));
    math::Vec f = sim.param().colvarres.K*w0*diff*derivatives[i];
    DEBUG(8, "Harmonic restr forces atom" << i+1 << ": " << math::v2s(f));
    (*atoms[i]).force(conf, topo, f);
    }
  }

  // entering harmonic linear restraining
  else if (sim.param().colvarres.colvarres == simulation::colvar_restr_linear_harmonic) {
    double shifting_value =  sim.param().distanceres.r_linear;
     for (int i=0; i<atoms.size(); i++) {
  
     }
  }
  return E;
}

int interaction::Colvar_Restraint_Interaction::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{
  // Contact number
  DEBUG(10, "Is contactnum_restrain empty: " << topo.contactnum_restraint().empty());
  if (!topo.contactnum_restraint().empty()) {
    DEBUG(8, "Entered Contact number restraining using COLVAR: ");
    for (int i=0; i< topo.contactnum_restraint().size(); i++) {
      interaction::Contactnum_Colvar *cv = new interaction::Contactnum_Colvar();
      (*cv).params = &topo.contactnum_restraint()[i];
      m_colvars.push_back(cv);
    }
  }

  // Distance restraint
  DEBUG(10, "Is distance_restraints_colvar empty: " << topo.distance_restraints_colvar().empty());
  if (!topo.distance_restraints_colvar().empty()) {
    DEBUG(8, "Entered Distance restraining using COLVAR: ");
    for (int i=0; i< topo.distance_restraints_colvar().size(); i++) {
      DEBUG(8, "distance restraint no: " << i+1 );
      interaction::Distance_Colvar *cv = new interaction::Distance_Colvar();
      (*cv).params = &topo.distance_restraints_colvar()[i];
      m_colvars.push_back(cv);
    }
  }
  
  // add more colvar types here

  
  for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
                                       to = m_colvars.end(); it != to; ++it) {
    if ((*it)->init(topo, conf, sim, os, quiet)) {
       os << "Colvar: error during initialisation of "
         << (*it)->name << "\n";

      io::messages.add("Error in colvar init",
                       "Colvar",
                       io::message::error);
    }
    Ctot0+=(*it)->targetvalue;
  }
    
  // if (!quiet) {
  //   os << "Colvar restraint interaction";
  //   os << std::endl;
  // }
  return 0;
}
