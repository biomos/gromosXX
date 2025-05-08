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
#include "../../interaction/special/colvar/perturbed_contactnum.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

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
  for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
    
    util::Virtual_Atom *firstatom = ((*it)->atoms)[0];
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
      
      (*it)->targetvalue=(1-l)*(*it)->targetvalueA + l*(*it)->targetvalueB;
      (*it)->w0=(1-l)*(*it)->w0A + l*(*it)->w0B;
        
    // calculate current value and derivatives
    if ((*it)->calculate(topo, conf, sim)) {
       std::cerr << "Colvar: error during calculation of "
         << (*it)->name << "\n";

      io::messages.add("Error in perturbed colvar calc",
                       "Colvar",
                       io::message::error);
    }
    Ctot+=(*it)->ct;
    conf.special().pertcolvarres.values.push_back((*it)->ct);
    
    // add forces
    if (!sum) {
      double E;
      E= apply_restraint(topo, conf, sim, (*it)->atoms, (*it)->derivatives, (*it)->ct, (*it)->targetvalue, (*it)->w0);
      conf.special().pertcolvarres.energies.push_back(E);
      Etot+=E;
      conf.current().energies.colvarres_energy[topo.atom_energy_group()
					    [(*firstatom).atom(0)]] += E;  
    }
    
    // get lambda derivatives
    double energy_derivativ, dpotdl;
    const double D_targetvalue = (*it)->targetvalueB - (*it)->targetvalueA;
    const double D_w0 = (*it)->w0B - (*it)->w0A;
    double diff = (*it)->ct - (*it)->targetvalue;
    
    if (sim.param().colvarres.colvarres == simulation::colvar_restr_harmonic) {
	  dpotdl = 0.5 * sim.param().colvarres.K * D_w0 * diff*diff
	  - sim.param().colvarres.K * (*it)->w0 * diff * D_targetvalue;
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
    for (std::vector<Colvar *>::iterator it = m_colvars.begin(),
        to = m_colvars.end(); it != to; ++it) {
      double E;
      // the individual weights are ignored when restraining to the sum
      E= apply_restraint(topo, conf, sim, (*it)->atoms, (*it)->derivatives, Ctot, Ctot0, 1.0);
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
                    math::VArray &derivatives, double &curr, double &target, double w0) {
                    
  double diff=target-curr;
  double E;
                  
  if (sim.param().colvarres.colvarres == simulation::colvar_restr_off) {
     
       std::cerr << "Colvar restraints are off, yet you ended up in colvar's apply_restraint .. something is wrong ..\n";

      io::messages.add("Error in colvar calc",
                       "Colvar",
                       io::message::error);
  }
  else if (sim.param().colvarres.colvarres == simulation::colvar_restr_harmonic) {
    E=0.5*sim.param().colvarres.K*w0*diff*diff;
//    if (sim.param().colvarres.write && ((sim.steps() - 1) % sim.param().colvarres.write) == 0)
//        std::cout << target << " - " << curr << " E: " << E << std::endl;
    
    for (int i=0; i<atoms.size(); i++) {
      math::Vec f = sim.param().colvarres.K*w0*diff*derivatives[i];
      (*atoms[i]).force(conf, topo, f);

        //if (sim.param().colvarres.virial) { 
        //  for (int a = 0; a < 3; ++a) {
	    //    for (int b = 0; b < 3; ++b) { 
	    //      conf.current().virial_tensor(a, b) += v(a) * f(b); 
	    //    }
        //  } 
        //}
      //if (sim.param().colvarres.write && ((sim.steps() - 1) % sim.param().colvarres.write) == 0) {
     //   if(math::abs(sim.param().colvarres.K*w0*diff*derivatives[i]) > 1) std::cout << "cvforces va " << i << " " << v2s(sim.param().colvarres.K*w0*diff*derivatives[i]) << ", w0 " << w0 <<std::endl;
     // }
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
  if (!topo.perturbed_contactnum_restraint().empty() and (sim.param().perturbation.perturbation)) {
    for (int i=0; i< topo.perturbed_contactnum_restraint().size(); i++) {
      interaction::Perturbed_Contactnum_Colvar *cv = new interaction::Perturbed_Contactnum_Colvar();
      (*cv).params = &topo.perturbed_contactnum_restraint()[i];
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
