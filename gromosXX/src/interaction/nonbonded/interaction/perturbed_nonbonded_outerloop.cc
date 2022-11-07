/**
 * @file perturbed_nonbonded_outerloop.cc
 * (template) methods of Perturbed_Nonbonded_Outerloop.
 */

#ifdef XXMPI
#include <mpi.h>
#endif

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_innerloop.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_innerloop.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h"

#include "../../../interaction/nonbonded/interaction_spec.h"

#include "../../../util/debug.h"
#include "../../../interaction/nonbonded/innerloop_template.h"

#include "../../interaction.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Perturbed_Nonbonded_Outerloop
::Perturbed_Nonbonded_Outerloop(Nonbonded_Parameter & nbp)
  : m_param(nbp)
{
}

//==================================================
// interaction loops
//==================================================

//==================================================
// the perturbed interaction (outer) loops
//==================================================

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_lj_crf_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Pairlist const & pairlist,
			     Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_perturbed_lj_crf_outerloop,
		       topo, conf, sim,
		       pairlist, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_lj_crf_outerloop(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      Pairlist const & pairlist,
			      Storage & storage)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);
  // Chris:
  // this is now (unfortunately) done in the innerloop, when we know which energy group we are in
  // innerloop.set_lambda(topo.lambda(), topo.lambda_exp());

  std::vector<unsigned int>::const_iterator j_it, j_to;
  unsigned int i = 0;
  unsigned int size_i = unsigned(topo.num_solute_atoms());
  // unsigned int size_i = unsigned(pairlist.size());

  DEBUG(6, "pert sr: " << size_i);

  for(i=0; i < size_i; ++i){
    
    for(j_it = pairlist[i].begin(),
	  j_to = pairlist[i].end();
	j_it != j_to;
	++j_it){
      
      DEBUG(10, "\tperturbed nonbonded_interaction: i "
	    << i << " j " << *j_it);
     // ANITA 
     // include sim such that we can decide if extendedTI is required 
     // innerloop.perturbed_lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
      innerloop.perturbed_lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity, sim);
    }
    
  }

  DEBUG(7, "end of function perturbed nonbonded interaction");  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_one_four_outerloop(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_perturbed_one_four_outerloop,
		       topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_one_four_outerloop(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				Storage & storage)
{
  DEBUG(7, "\tcalculate perturbed 1,4-interactions");
  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);

  // Chris:
  // this is now (unfortunately) done in the innerloop, when we know which energy group we are in
  // innerloop.set_lambda(topo.lambda(), topo.lambda_exp());
  
  topology::excl_cont_t::value_type::const_iterator it, to;
  std::map<unsigned int, topology::Perturbed_Atom>::const_iterator 
    mit=topo.perturbed_solute().atoms().begin(), 
    mto=topo.perturbed_solute().atoms().end();
  
  for(; mit!=mto; ++mit){
    it = mit->second.one_four_pair().begin();
    to = mit->second.one_four_pair().end();
    
    for( ; it != to; ++it){

      innerloop.perturbed_one_four_interaction_innerloop
	(topo, conf, mit->second.sequence_number(), *it, periodicity, sim); // ANITA: include sim

    } // loop over 1,4 pairs
  } // loop over solute atoms
}

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_RF_excluded_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_perturbed_RF_excluded_outerloop,
		       topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_RF_excluded_outerloop(topology::Topology & topo,
				   configuration::Configuration & conf,
				   simulation::Simulation & sim,
				   Storage & storage)
{

  DEBUG(7, "\tcalculate perturbed excluded RF interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);
  // Chris:
  // this is now (unfortunately) done in the innerloop, when we know which energy group we are in
  // innerloop.set_lambda(topo.lambda(), topo.lambda_exp());

  std::map<unsigned int, topology::Perturbed_Atom>::const_iterator
    mit=topo.perturbed_solute().atoms().begin(),
    mto=topo.perturbed_solute().atoms().end();

  DEBUG(9, "\tSize of perturbed atoms " 
	<< unsigned(topo.perturbed_solute().atoms().size()));
  
  for(; mit!=mto; ++mit){
    innerloop.perturbed_RF_excluded_interaction_innerloop(topo, conf, mit, periodicity, sim); //ANITA
  }
}

void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_electric_field_outerloop(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim, 
		   PairlistContainer const & pairlist,
                   PairlistContainer const & perturbed_pairlist,
                   Storage & storage,
                   Storage & storage_lr,
                   int rank)
{
  SPLIT_PERT_INNERLOOP(_perturbed_electric_field_outerloop, topo, conf, sim, 
                  pairlist, perturbed_pairlist, storage, storage_lr, rank);
}
/**
 * helper function to calculate polarisation, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_electric_field_outerloop(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
		    PairlistContainer const & pairlist,
                    PairlistContainer const & perturbed_pairlist,
		    Storage & storage,
                    Storage & storage_lr,
                    int rank)
{  
  DEBUG(7, "\tcalculate polarisation (electric field outerloop)");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  // unperturbed innerloop
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  // perturbed innerloop
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> p_innerloop(m_param);
  p_innerloop.init(sim);
  // Chris:
  // this is now (unfortunately) done in the innerloop, when we know which energy group we are in
  // p_innerloop.set_lambda(topo.lambda(), topo.lambda_exp());
  
  unsigned int i = 0;
  unsigned int size_i = unsigned(pairlist.size());
  unsigned int size_lr = size_i;
  DEBUG(11, "outerloop pairlist size " << size_i);

  unsigned int end = size_i;
  unsigned int end_lr = size_lr;

  math::VArray e_el_new(topo.num_atoms());
  
#ifdef XXMPI
  // because we need some place to reduce the field to
  math::VArray e_el_master(topo.num_atoms());
#endif

  double minfield = sim.param().polarise.minfield;
  const double minfield_param = minfield;
  double maxfield = 0.0;
  int turni = 0;

#ifdef XXMPI
  // broadcast posV to slaves. We only have to do this here at the very first step because
  // posV is also broadcasted at the end of every electric field iteration.
  if (sim.mpi && sim.steps() == 0) {
    MPI_Bcast(&conf.current().posV(0)(0), conf.current().posV.size() * 3, MPI::DOUBLE, sim.mpiControl().masterID, sim.mpiControl().comm);
  }
#endif

  // longrange ?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){

    // loop over all molecules in longrange pairlist
    for(i=0; i < end_lr; ++i){
      // perturbed solute longrange
      for(std::vector<unsigned int>::const_iterator 
              j_it = perturbed_pairlist.solute_long[i].begin(),
              j_to = perturbed_pairlist.solute_long[i].end();
              j_it != j_to; ++j_it){
        
        math::Vec e_eli_lr, e_elj_lr;
        
        p_innerloop.perturbed_electric_field_innerloop(topo, conf, i, *j_it,
                e_eli_lr, e_elj_lr, periodicity);
        
        storage_lr.electric_field[i] += e_eli_lr;
        storage_lr.electric_field[*j_it] += e_elj_lr;
      }
      // solute longrange
      for(std::vector<unsigned int>::const_iterator 
              j_it = pairlist.solute_long[i].begin(),
              j_to = pairlist.solute_long[i].end();
              j_it != j_to; ++j_it){
        
        math::Vec e_eli_lr, e_elj_lr;

        innerloop.electric_field_innerloop(topo, conf, i, *j_it,
                e_eli_lr, e_elj_lr, periodicity);

        storage_lr.electric_field[i] += e_eli_lr;
        storage_lr.electric_field[*j_it] += e_elj_lr;
      }
      // solvent longrange
      for(std::vector<unsigned int>::const_iterator 
              j_it = pairlist.solvent_long[i].begin(),
              j_to = pairlist.solvent_long[i].end();
              j_it != j_to; ++j_it){
        
        math::Vec e_eli_lr, e_elj_lr;

        innerloop.electric_field_innerloop(topo, conf, i, *j_it,
                e_eli_lr, e_elj_lr, periodicity);

        storage_lr.electric_field[i] += e_eli_lr;
        storage_lr.electric_field[*j_it] += e_elj_lr;
      }
    }
#ifdef XXMPI
    if (sim.mpi) {
      // reduce the longrange electric field to some temp. variable and then set this
      // variable to the longrange electric field on the master. The lr e field
      // is only needed on the master node
      if (rank) {
        MPI_Reduce(&storage_lr.electric_field(0)(0), NULL,
                             storage_lr.electric_field.size() * 3, MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
      } else {
        MPI_Reduce(&storage_lr.electric_field(0)(0), &e_el_master(0)(0),
                             storage_lr.electric_field.size() * 3, MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
        storage_lr.electric_field = e_el_master;
      }
    }
#endif
  }

  // shortrange
  while (minfield >= minfield_param) {

    maxfield = 0.0;
    e_el_new = 0.0;
    
#ifdef XXMPI
    // again set the temporary variable to 0 as we need it again for 
    // the short range eletric field
    if (sim.mpi)
      e_el_master = 0.0;
#endif
    
    // loop over all molecules in shortrange pairlist
    for(i=0; i < end; ++i){
       // perturbed solute short
      for(std::vector<unsigned int>::const_iterator 
              j_it = perturbed_pairlist.solute_short[i].begin(),
	      j_to = perturbed_pairlist.solute_short[i].end();
              j_it != j_to; ++j_it){

        math::Vec e_eli, e_elj;

        p_innerloop.perturbed_electric_field_innerloop(topo, conf, 
           i, *j_it, e_eli, e_elj, periodicity);

        e_el_new(i) += e_eli;
        e_el_new(*j_it) += e_elj;
      }
      // solute short
      for(std::vector<unsigned int>::const_iterator 
              j_it = pairlist.solute_short[i].begin(),
	      j_to = pairlist.solute_short[i].end();
              j_it != j_to; ++j_it){

        math::Vec e_eli, e_elj;

        innerloop.electric_field_innerloop(topo, conf, 
           i, *j_it, e_eli, e_elj, periodicity);

        e_el_new(i) += e_eli;
        e_el_new(*j_it) += e_elj;
      }
      // solvent short
      for(std::vector<unsigned int>::const_iterator 
              j_it = pairlist.solvent_short[i].begin(),
	      j_to = pairlist.solvent_short[i].end();
              j_it != j_to; ++j_it){

        math::Vec e_eli, e_elj;

        innerloop.electric_field_innerloop(topo, conf, 
           i, *j_it, e_eli, e_elj, periodicity);

        e_el_new(i) += e_eli;
        e_el_new(*j_it) += e_elj;
      }
    }

#ifdef XXMPI
    // also reduce the shortrange electric field the same way as the longrange
    // electric field
    if (sim.mpi) {
      if (rank) {
        MPI_Reduce(&e_el_new(0)(0), NULL, e_el_new.size() * 3, MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
      } else {
        MPI_Reduce(&e_el_new(0)(0), &e_el_master(0)(0), e_el_new.size() * 3, MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
        e_el_new = e_el_master;
      }
    }
#endif

    if (rank == 0) {
      
      for (i=0; i<topo.num_atoms(); ++i) {
        double alpha = topo.polarisability(i);
        double damp_lev =  topo.damping_level(i);
        const double damp_pow = topo.damping_power(i);
        
        // (1-l)^n * X + l^n * X
        if (topo.is_perturbed(i)) {
	  const double lambda = topo.individual_lambda(simulation::crf_lambda)
	    [topo.atom_energy_group()[i]][topo.atom_energy_group()[i]];
	  
          alpha = pow((1-lambda), topo.lambda_exp()) * 
	    topo.perturbed_solute().atoms()[i].A_polarisability()
	    + pow(lambda, topo.lambda_exp()) * 
	    topo.perturbed_solute().atoms()[i].B_polarisability();
          damp_lev = pow((1-lambda), topo.lambda_exp()) * 
	    topo.perturbed_solute().atoms()[i].A_damping_level()
	    + pow(lambda, topo.lambda_exp()) * 
	    topo.perturbed_solute().atoms()[i].B_damping_level();
        }
        if(topo.is_polarisable(i)){
          e_el_new(i) += storage_lr.electric_field(i);
          
          //delta r
          math::Vec delta_r;

          //////////////////////////////////////////////////
          // implementation of polarisability damping
          /////////////////////////////////////////////////

          if (sim.param().polarise.damp) { // damp the polarisability
            const double e_i = sqrt(math::abs2(e_el_new(i))),
                    e_0 = damp_lev;
            if (e_i <= e_0)
              delta_r = (alpha / topo.coscharge(i)) * e_el_new(i);
            else {
              const double p = damp_pow;
              delta_r = alpha * e_0 *
                        (p + 1.0 - pow(e_0/e_i, p)) /
		        (p*topo.coscharge(i) * e_i) * e_el_new(i);
            }
          } else { // no damping
            delta_r = (alpha / topo.coscharge(i)) * e_el_new(i);
          }
          // store the new position
          conf.current().posV(i) = delta_r;

          // calculation of convergence criterium
          for(int j=0; j<3; ++j) {
            double delta_e = fabs(storage.electric_field(i)(j)-e_el_new(i)(j));
            if (delta_e > maxfield) {
              maxfield = delta_e;
            }
          }
        }
        storage.electric_field(i) = e_el_new(i);
      }
    } // end master only    
    
    turni++;
    minfield = maxfield;

#ifdef XXMPI
    // broadcast the new posV and also the convergence criterium (minfield)
    // to the slaves. Otherwise they don't know when to stop.
    if (sim.mpi) {
      MPI_Bcast(&conf.current().posV(0)(0), conf.current().posV.size() * 3, MPI::DOUBLE, sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Bcast(&minfield, 1, MPI::DOUBLE, sim.mpiControl().masterID, sim.mpiControl().comm);
    }
#endif
    DEBUG(11, "\trank: " << rank << " minfield: "<<minfield<<" iteration round: "<<turni);
  }
}
/**
 * helper function to calculate self energy, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_self_energy_outerloop(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim, 
                   Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_perturbed_self_energy_outerloop, topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_self_energy_outerloop(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
		    Storage & storage)
{  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop_p(m_param);
  innerloop_p.init(sim);
  // Chris:
  // this is now (unfortunately) done in the innerloop, when we know which energy group we are in
  // innerloop_p.set_lambda(topo.lambda(), topo.lambda_exp());

  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    if(topo.is_polarisable(i)){
      DEBUG (10, "\tperturbed self energy outerloop");
      if (topo.is_perturbed(i)) {
        innerloop_p.perturbed_self_energy_innerloop(topo, conf, i, 
                                        storage, periodicity);
      }
      else {
        innerloop.self_energy_innerloop(topo, conf, i, 
                                        storage, periodicity);
      }
    }
  }
}


