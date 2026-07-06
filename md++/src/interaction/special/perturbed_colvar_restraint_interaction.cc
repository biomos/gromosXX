/**
 * @file perturbed_colvar_restraint_interaction.cc
 * template methods of Perturbed_Colvar_Restraint_Interaction
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

#include "../../interaction/special/perturbed_colvar_restraint_interaction.h"
#include "../../interaction/special/colvar/colvar.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

namespace {

bool colvar_type_matches(const std::string &stored, const std::string &requested) {
  if (stored == requested) return true;

  // Accepted aliases for perturbed coordination-number restraints.
  if (requested == "PERTCOORDNUM") {
    return stored == "PERTCOORDNUM" ||
           stored == "PERT_COORDNUM" ||
           stored == "PERTURBED_COORDNUM" ||
           stored == "PERTURBED_COORDNUM";
  }

  return false;
}

unsigned int count_colvar_specs(const simulation::Simulation &sim,
                                const std::string &type) {
  unsigned int count = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = sim.param().colvarres.bias_specs.begin(),
       to = sim.param().colvarres.bias_specs.end(); it != to; ++it) {
    if (colvar_type_matches(it->type, type)) ++count;
  }
  return count;
}

double colvar_k_for_type_instance(const simulation::Simulation &sim,
                                  const std::string &type,
                                  const unsigned int instance) {
  unsigned int count = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = sim.param().colvarres.bias_specs.begin(),
       to = sim.param().colvarres.bias_specs.end(); it != to; ++it) {
    if (!colvar_type_matches(it->type, type)) continue;
    if (count == instance) return it->k;
    ++count;
  }

  io::messages.add("COLVARRES: missing bias specification for perturbed colvar type "
                   + type + ". Force constant set to 0.0.",
                   "Perturbed_Colvar_Restraint_Interaction",
                   io::message::error);
  return 0.0;
}

} // anonymous namespace

interaction::Perturbed_Colvar_Restraint_Interaction::~Perturbed_Colvar_Restraint_Interaction(){
  for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
        delete *it;
  }
}


/**
 * calculate colvar restraint interactions
 */
int interaction::Perturbed_Colvar_Restraint_Interaction
::calculate_interactions(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  m_timer.start(sim);
  double Ctot=0;
  double Etot=0;
  conf.special().pertcolvarres.energies.clear();
  conf.special().pertcolvarres.values.clear();
  
  // loop over all colvars 
  for (size_t cv_index = 0; cv_index < m_colvars.size(); ++cv_index) {
    Colvar *cv = m_colvars[cv_index];
    const double K = m_force_constants[cv_index];
    
    util::Virtual_Atom *firstatom = (cv->atoms)[0];
    // the first atom of the first virtual atom 
    // determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double l = topo.individual_lambda(simulation::colvarres_lambda)
      [topo.atom_energy_group()[(*firstatom).atom(0)]]
      [topo.atom_energy_group()[(*firstatom).atom(0)]];
    const double l_deriv = topo.individual_lambda_derivative
      (simulation::colvarres_lambda)
      [topo.atom_energy_group()[(*firstatom).atom(0)]]
      [topo.atom_energy_group()[(*firstatom).atom(0)]];
      
    cv->targetvalue=(1-l)*cv->targetvalueA + l*cv->targetvalueB;
    cv->w0=(1-l)*cv->w0A + l*cv->w0B;
        
    // calculate current value and derivatives
    if (cv->calculate(topo, conf, sim)) {
       std::cerr << "Colvar: error during calculation of "
         << cv->name << "\n";

      io::messages.add("Error in perturbed colvar calc",
                       "Colvar",
                       io::message::error);
    }
    Ctot+=cv->ct;
    conf.special().pertcolvarres.values.push_back(cv->ct);
    
    // add forces
    if (!sum) {
      double E;
      E= apply_restraint(topo, conf, sim, cv->atoms, cv->derivatives,
                         cv->ct, cv->targetvalue, cv->w0, K);
      conf.special().pertcolvarres.energies.push_back(E);
      Etot+=E;
      conf.current().energies.colvarres_energy[topo.atom_energy_group()
                                            [(*firstatom).atom(0)]] += E;  
    }
    
    // get lambda derivatives
    double energy_derivativ = 0.0, dpotdl = 0.0;
    const double D_targetvalue = cv->targetvalueB - cv->targetvalueA;
    const double D_w0 = cv->w0B - cv->w0A;
    double diff = cv->ct - cv->targetvalue;
    
    if (sim.param().colvarres.colvarres == simulation::colvar_restr_harmonic) {
          dpotdl = 0.5 * K * D_w0 * diff*diff
          - K * cv->w0 * diff * D_targetvalue;
        }
        energy_derivativ = l_deriv * dpotdl;
    DEBUG(10, "PERTCOLVARRES Energy derivative : " << energy_derivativ);

    conf.current().perturbed_energy_derivatives.
      colvarres_energy[topo.atom_energy_group()
                                            [(*firstatom).atom(0)]] += 
      energy_derivativ;
  }
  
  conf.special().pertcolvarres.totv=Ctot;
  
  // not implemented yet, restrain to the sum of all colvar restraints
  if (sum) {
    const double K_sum = m_force_constants.empty() ? 0.0 : m_force_constants[0];
    for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
      double E;
      // the individual weights are ignored when restraining to the sum
      E= apply_restraint(topo, conf, sim, (*it)->atoms, (*it)->derivatives,
                         Ctot, Ctot0, 1.0, K_sum);
      conf.special().pertcolvarres.energies.push_back(E);
      Etot+=E;
      // add to the energy group of the first atom of first list
      util::Virtual_Atom *firstatom = ((*it)->atoms)[0];
      conf.current().energies.colvarres_energy[topo.atom_energy_group()
                                            [(*firstatom).atom(0)]] += E;
    }
  }
  conf.special().pertcolvarres.tote=Etot;
  m_timer.stop();    
  return 0;
}


// apply forces and calculate potential using the biasing function specified in the COLVARRES block
double interaction::Perturbed_Colvar_Restraint_Interaction
  ::apply_restraint(topology::Topology & topo,
 configuration::Configuration & conf,simulation::Simulation & sim, std::vector< util::Virtual_Atom* > atoms, 
                    math::VArray &derivatives, double &curr, double &target, double w0,
                    double K) {
                    
  double diff=target-curr;
  double E = 0.0;
                  
  if (sim.param().colvarres.colvarres == simulation::colvar_restr_off) {
     
       std::cerr << "Colvar restraints are off, yet you ended up in colvar's apply_restraint .. something is wrong ..\n";

      io::messages.add("Error in colvar calc",
                       "Colvar",
                       io::message::error);
  }
  else if (sim.param().colvarres.colvarres == simulation::colvar_restr_harmonic) {
    E=0.5*K*w0*diff*diff;
    
    for (size_t i=0; i<atoms.size(); i++) {
      math::Vec f = K*w0*diff*derivatives[i];
      (*atoms[i]).force(conf, topo, f);
    }
  }
  return E;
}

int interaction::Perturbed_Colvar_Restraint_Interaction::init(topology::Topology &topo, 
                     configuration::Configuration &conf,
                     simulation::Simulation &sim,
                     std::ostream &os,
                     bool quiet) 
{
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