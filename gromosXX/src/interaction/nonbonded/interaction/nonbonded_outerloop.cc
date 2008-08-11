/**
 * @file nonbonded_outerloop.cc
 * template methods of Nonbonded_Outerloop.
 */

#ifdef XXMPI
#include <mpi.h>
#endif

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/latticesum.h>
#include <configuration/mesh.h>
#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>

#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>

#include <interaction/nonbonded/interaction_spec.h>

#include <util/debug.h>
#include <interaction/nonbonded/innerloop_template.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Nonbonded_Outerloop
::Nonbonded_Outerloop(Nonbonded_Parameter &nbp)
  : m_param(nbp)
{
}

//==================================================
// interaction loops
//==================================================

void interaction::Nonbonded_Outerloop
::lj_crf_outerloop(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim, 
                   Pairlist const & pairlist_solute,
                   Pairlist const & pairlist_solvent,
		   Storage & storage)
{
  SPLIT_INNERLOOP(_lj_crf_outerloop, topo, conf, sim,
                  pairlist_solute, pairlist_solvent, storage);
}


/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_lj_crf_outerloop(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
	            Pairlist const & pairlist_solute,
                    Pairlist const & pairlist_solvent,
		    Storage & storage)
{  
  DEBUG(7, "\tcalculate interactions");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  /*
    variables for a OMP parallelizable loop.
    outer index has to be integer...
  */
  std::vector<unsigned int>::const_iterator j_it, j_to;

  unsigned int size_i = unsigned(pairlist_solute.size());
  DEBUG(10, "outerloop pairlist size " << size_i);
  
  const unsigned int end = topo.num_solute_atoms();
  
  unsigned int i;
  for(i=0; i < end; ++i){
    for(j_it = pairlist_solute[i].begin(),
	  j_to = pairlist_solute[i].end();
	j_it != j_to;
	++j_it){
      
      DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << *j_it);

      // shortrange, therefore store in simulation.system()
      innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }
  
  if (sim.param().force.special_loop == simulation::special_loop_spc) { // special solvent loop
    // solvent - solvent with spc innerloop...
    for(; i < size_i; i += 3){ // use every third pairlist (OW's)
      for(j_it = pairlist_solvent[i].begin(),
              j_to = pairlist_solvent[i].end();
      j_it != j_to;
      j_it += 3){ // use every third atom (OW) in pairlist i
        
        DEBUG(10, "\tspc_nonbonded_interaction: i " << i << " j " << *j_it);
        
        innerloop.spc_innerloop(topo, conf, i, *j_it, storage, periodicity);
      }
    }
  } else { // normal solvent loop
    for(; i < size_i; ++i){
      for(j_it = pairlist_solvent[i].begin(),
              j_to = pairlist_solvent[i].end();
      j_it != j_to;
      ++j_it){
        
        DEBUG(10, "\tsolvent_nonbonded_interaction: i " << i << " j " << *j_it);
        
        innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
      }
    }
  }
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
void interaction::Nonbonded_Outerloop
::one_four_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage)
{
  SPLIT_INNERLOOP(_one_four_outerloop, topo, conf, sim, storage);
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_one_four_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage)
{
  DEBUG(7, "\tcalculate 1,4-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  
  std::set<int>::const_iterator it, to;
  unsigned int const num_solute_atoms = topo.num_solute_atoms();
  unsigned int i;
  
  for(i=0; i < num_solute_atoms; ++i){
    it = topo.one_four_pair(i).begin();
    to = topo.one_four_pair(i).end();
    
    for( ; it != to; ++it){

      innerloop.one_four_interaction_innerloop(topo, conf, i, *it, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
void interaction::Nonbonded_Outerloop
::cg_exclusions_outerloop(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  Storage & storage)
{
  SPLIT_INNERLOOP(_cg_exclusions_outerloop, topo, conf, sim, storage);
}

template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_cg_exclusions_outerloop(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  Storage & storage)
{
  DEBUG(7, "\tcalculate 1,4-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  
  std::set<int>::const_iterator cg2_it, cg2_to;

  for(unsigned int cg1=0; cg1<topo.num_solute_chargegroups(); ++cg1){
    
    cg2_it = topo.chargegroup_exclusion(cg1).begin();
    cg2_to = topo.chargegroup_exclusion(cg1).end();
    
    for( ; cg2_it != cg2_to; ++cg2_it){

      // only once...
      if (cg1 > (unsigned) *cg2_it) continue;
      
      for(int a1 = topo.chargegroup(cg1); a1 < topo.chargegroup(cg1 + 1); ++a1){
	for(int a2 = topo.chargegroup(*cg2_it); a2 < topo.chargegroup(*cg2_it + 1); ++a2){
	  
	  // std::cout << "cg1=" << cg1 << " cg2=" << *cg2_it 
	  // << " a1=" << a1 << " a2=" << a2 << std::endl;
	  
	  if (a1 >= a2) continue;
	  
	  if (topo.exclusion(a1).find(a2) != topo.exclusion(a1).end()) continue;
	  
	  if (topo.one_four_pair(a1).find(a2) != topo.one_four_pair(a1).end()){
	    // std::cout << "\t1,4" << std::endl;
	    innerloop.one_four_interaction_innerloop(topo, conf, a1, a2, periodicity);
	  }
	  else{
	    // std::cout << "\tstandard interaction" << std::endl;
	    innerloop.lj_crf_innerloop(topo, conf, a1, a2, storage, periodicity);
	  }
	} // atoms of cg 2
      } // atoms of cg 1

    } // cg 2 (excluded from cg 1)
  } // solute cg's
}  


/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
void interaction::Nonbonded_Outerloop
::RF_excluded_outerloop(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			Storage & storage)
{
  /*
  if (sim.param().force.interaction_function !=
      simulation::lj_crf_func &&
      sim.param().force.interaction_function !=
      simulation::pol_lj_crf_func){
    io::messages.add("Nonbonded_Outerloop",
		     "RF excluded term for non-lj_crf_func called",
		     io::message::error);
  }
   */
  
  SPLIT_INNERLOOP(_RF_excluded_outerloop, topo, conf, sim, storage);  
}

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_RF_excluded_outerloop(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			Storage & storage)
{
  
  DEBUG(7, "\tcalculate RF excluded interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  
  int i, num_solute_atoms = topo.num_solute_atoms();
  
  for(i=0; i<num_solute_atoms; ++i){
    innerloop.RF_excluded_interaction_innerloop(topo, conf, i, periodicity);
  } // loop over solute atoms

  // Solvent
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();
  cg_it += topo.num_solute_chargegroups();

  for(; cg_it != cg_to; ++cg_it){
    innerloop.RF_solvent_interaction_innerloop(topo, conf, cg_it, periodicity);
  } // loop over solvent charge groups
}  

void interaction::Nonbonded_Outerloop
::self_energy_outerloop(topology::Topology & topo,
		        configuration::Configuration & conf,
		        simulation::Simulation & sim, 
			Storage & storage)
{
  SPLIT_INNERLOOP(_self_energy_outerloop, topo, conf, sim, storage);
}

/**
 * helper function to calculate self energy, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_self_energy_outerloop(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
                    Storage & storage)
{  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    if(topo.is_polarizable(i)){
    
      DEBUG(10, "\tself energy: i " << i);
      innerloop.self_energy_innerloop (
           topo, conf, i, storage, periodicity);
    }
  }
}


void interaction::Nonbonded_Outerloop
::electric_field_outerloop(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim, 
		   PairlistContainer const & pairlist,
                   Storage & storage,
                   Storage & storage_lr, int rank)
{
  SPLIT_INNERLOOP(_electric_field_outerloop, topo, conf, sim, 
                  pairlist, storage, storage_lr, rank);
}
/**
 * helper function to calculate polarization, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_electric_field_outerloop(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
		    PairlistContainer const & pairlist,
		    Storage & storage,
                    Storage & storage_lr,
                    int rank)
{  
  DEBUG(7, "\tcalculate polarization (electric field outerloop)");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  unsigned int i;
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

  double minfield = sim.param().polarize.minfield;
  const double minfield_param = minfield;
  double maxfield;
  int turni = 0;

#ifdef XXMPI
  // broadcast posV to slaves. We only have to do this here at the very first step because
  // posV is also broadcasted at the end of every electric field iteration.
  if (sim.mpi && sim.steps() == 0) {
    MPI::COMM_WORLD.Bcast(&conf.current().posV(0)(0), conf.current().posV.size() * 3, MPI::DOUBLE, 0);
  }
#endif
 
  // longrange ?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){

    // loop over all molecules in longrange pairlist
    for(i=0; i < end_lr; ++i){
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
        MPI::COMM_WORLD.Reduce(&storage_lr.electric_field(0)(0), NULL,
                             storage_lr.electric_field.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
      } else {
        MPI::COMM_WORLD.Reduce(&storage_lr.electric_field(0)(0), &e_el_master(0)(0),
                             storage_lr.electric_field.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
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
        MPI::COMM_WORLD.Reduce(&e_el_new(0)(0), NULL, e_el_new.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
      } else {
        MPI::COMM_WORLD.Reduce(&e_el_new(0)(0), &e_el_master(0)(0), e_el_new.size() * 3, MPI::DOUBLE, MPI::SUM, 0);
        e_el_new = e_el_master;
      }
    }
#endif

    if (rank == 0) {
      for (i=0; i<topo.num_atoms(); ++i) {
        if(topo.is_polarizable(i)){
          e_el_new(i) += storage_lr.electric_field(i);

          //delta r
          math::Vec delta_r;
        
          //////////////////////////////////////////////////
          // implementation of polarizability damping
          /////////////////////////////////////////////////
        
          if (sim.param().polarize.damp) { // damp the polarizability
            const double e_i = sqrt(math::abs2(e_el_new(i))),
                         e_0 = topo.damping_level(i);
            if (e_i <= e_0) 
              delta_r = (topo.polarizability(i) / topo.coscharge(i)) * e_el_new(i);
            else {
              const double p = topo.damping_power(i);
              delta_r = topo.polarizability(i) * e_0 / p * 
                        (p + 1.0 - pow(e_0/e_i, p)) / 
                        (topo.coscharge(i) * e_i) * e_el_new(i);
            }
          } else { // no damping
            delta_r = (topo.polarizability(i) / topo.coscharge(i)) * e_el_new(i);
          }
          // store the new position
          conf.current().posV(i) = delta_r;
        
          // calculation of convergence criterium
          for(int j=0; j<3; ++j) {
            double delta_e = fabs(storage.electric_field(i)(j)-e_el_new(i)(j))* 7.911492226513023 * 0.1;
            if (delta_e > maxfield) {
              maxfield = delta_e;
            }
          }
        }
        storage.electric_field(i) = e_el_new(i);
      }
    }
    turni++;
    minfield = maxfield;

#ifdef XXMPI
    // broadcast the new posV and also the convergence criterium (minfield)
    // to the slaves. Otherwise they don't know when to stop.
    if (sim.mpi) {
      MPI::COMM_WORLD.Bcast(&conf.current().posV(0)(0), conf.current().posV.size() * 3, MPI::DOUBLE, 0);
      MPI::COMM_WORLD.Bcast(&minfield, 1, MPI::DOUBLE, 0);
    }
#endif
    DEBUG(11, "\trank: " << rank << " minfield: "<<minfield<<" iteration round: "<<turni);
  }
}

void interaction::Nonbonded_Outerloop
::ls_real_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Pairlist const & pairlist_solute,
        Pairlist const & pairlist_solvent,
        Storage & storage, int rank, int size) {
  SPLIT_INNERLOOP(_ls_real_outerloop, topo, conf, sim,
          pairlist_solute, pairlist_solvent, storage, rank, size);
}

template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_ls_real_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Pairlist const & pairlist_solute,
        Pairlist const & pairlist_solvent,
        Storage & storage, int rank, int size)
{  
  
  DEBUG(7, "\tcalculate LS real space interactions");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  /*
    variables for a OMP parallelizable loop.
    outer index has to be integer...
  */
  std::vector<unsigned int>::const_iterator j_it, j_to;
  std::set<int>::const_iterator ex_it, ex_to;

  unsigned int size_i = unsigned(pairlist_solute.size());
  DEBUG(10, "outerloop pairlist size " << size_i);
  
  unsigned int end = topo.num_solute_atoms();
    
  unsigned int i;
  for(i=0; i < end; i++){
    for(j_it = pairlist_solute[i].begin(),
	  j_to = pairlist_solute[i].end();
	j_it != j_to;
	++j_it){
      
      DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << *j_it);

      // shortrange, therefore store in simulation.system()
      innerloop.lj_ls_real_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  } 
  
  for(; i < size_i; i++){
    for(j_it = pairlist_solvent[i].begin(),
            j_to = pairlist_solvent[i].end();
    j_it != j_to;
    ++j_it){
      
      DEBUG(10, "\tsolvent_nonbonded_interaction: i " << i << " j " << *j_it);
      
      innerloop.lj_ls_real_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }
  
  // loop over all exclusions as they have reduced LS interactions in real space
  // parrallelization using stride. Then MPI should work
  DEBUG(9, "U_eta due to excluded solute pairs...");
  const unsigned int size_int = topo.num_solute_atoms();
  for(unsigned int i = rank; i < size_int; i+=size) {
    for(ex_it = topo.exclusion(i).begin(),
        ex_to = topo.exclusion(i).end();
    ex_it != ex_to;
    ++ex_it){
      
      DEBUG(10, "\texcluded_nonbonded_interaction: i " << i << " j " << *ex_it);
      innerloop.ls_real_excluded_innerloop(topo, conf, i, *ex_it, storage, periodicity);
    }
  }
  DEBUG(9, "U_eta due to excluded solvent pairs...");
  const unsigned int num_cg = topo.num_chargegroups();
  for(unsigned int i = topo.num_solute_chargegroups() + rank; i < num_cg; i+=size) {
    for (int a1 = topo.chargegroup(i),
            a_to = topo.chargegroup(i + 1);
            a1 < a_to; ++a1) {
      for (int a2 = a1 + 1; a2 < a_to; ++a2) {
        DEBUG(10, "\texcluded_nonbonded_interaction: i " << a1 << " j " << a2);
        innerloop.ls_real_excluded_innerloop(topo, conf, a1, a2, storage, periodicity);
      }
    }
  }
   
}

void interaction::Nonbonded_Outerloop
::ls_ewald_kspace_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size) {
  SPLIT_INNERLOOP(_ls_ewald_kspace_outerloop, topo, conf, sim,
          storage, rank, size);
}
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_ls_ewald_kspace_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size)
{  
  DEBUG(7, "\tcalculate interactions in k-space (Ewald)");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  const double volume = math::volume(conf.current().box, t_interaction_spec::boundary_type);
  const double eps_volume_i = math::eps0_i / volume;
  
  DEBUG(10, "\t\teps_volume_i: " << eps_volume_i);
  const std::vector<configuration::KSpace_Element> & kspace = conf.lattice_sum().kspace;
  std::vector<configuration::KSpace_Element>::const_iterator it = kspace.begin(),
          to = kspace.end();
  
  const unsigned int num_atoms = topo.num_atoms();
  // a copy of the current positions because of gathering
  math::VArray r = conf.current().pos;
  
  // storage for the k space energy
  double energy = 0.0;
  // and force
  math::VArray f(num_atoms);
  f = 0.0;
  
  // Do we have to gather here ?
  for(unsigned int i = 0; i < num_atoms; ++i) {
    periodicity.put_into_positive_box(r(i));
    DEBUG(11, "r(" << i <<") in box: " << math::v2s(r(i)));
  }
  
  // cache for sin(k . r) and cos(k .  r) terms
  math::SArray sin_kr(num_atoms, 0.0), cos_kr(num_atoms, 0.0);
  
  // on fly we can calculate the methodology dependent A2 term
  double a2_tilde = 0.0;
  
  // loop over k space
  for(; it != to; ++it) {
    DEBUG(12, "k: " << math::v2s(it->k));
    double C_k = 0.0;
    double S_k = 0.0;
    // loop over atoms, calculate C/S_k and cache the sins and coses
    for(unsigned int i = 0; i < num_atoms; ++i) {
      const double r_dot_k = math::dot(it->k, r(i));
      sin_kr(i) = sin(r_dot_k);
      cos_kr(i) = cos(r_dot_k);
      C_k += topo.charge(i) * cos_kr(i);
      S_k += topo.charge(i) * sin_kr(i);
    }
    
    // calculate the force component from k
    for(unsigned int i = 0; i < num_atoms; ++i) {
      f(i) += it->k * (it->k2i_gammahat * (C_k * sin_kr(i) - S_k * cos_kr(i)));
    }
    // and the energy component from k
    energy += it->k2i_gammahat * (C_k * C_k + S_k * S_k);
    // add to the A2 term
    a2_tilde += it->k2i_gammahat;
  }
 
  // loop again over atoms and store force
  for(unsigned int i = 0; i < num_atoms; ++i) {
    // calculate the force
    f(i) *= topo.charge(i) * eps_volume_i;
    DEBUG(12, "force f(" << i << "): " << math::v2s(f(i)));
    // add the force
    storage.force(i) += f(i);
  }
  
  // scale the energy
  energy *= eps_volume_i * 0.5; 
  DEBUG(8, "Ewald k-space energy: " << energy);
  storage.energies.ls_kspace_total = energy;
  
  // scale the a2 term
  a2_tilde *= 4.0 * math::Pi / volume;
  conf.lattice_sum().a2_tilde = a2_tilde;
}

void interaction::Nonbonded_Outerloop
::ls_p3m_kspace_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size,
        util::Algorithm_Timer & timer) {
  SPLIT_INNERLOOP(_ls_p3m_kspace_outerloop, topo, conf, sim,
          storage, rank, size, timer);
}
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_ls_p3m_kspace_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size,
        util::Algorithm_Timer & timer)
{  
  DEBUG(7, "\tcalculate interactions in k-space (P3M)");  
    

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  math::VArray r = conf.current().pos;
  // Do we have to gather here ?
  const unsigned int num_atoms = topo.num_atoms();
  for (unsigned int i = 0; i < num_atoms; ++i) {
    periodicity.put_into_positive_box(r(i));
    DEBUG(11, "r(" << i << ") in box: " << math::v2s(r(i)));
  }
  
  // decompose into domains
  if (sim.mpi)
    interaction::Lattice_Sum::decompose_into_domains<configuration::ParallelMesh>(topo, conf, sim, r, size);
  else
    interaction::Lattice_Sum::decompose_into_domains<configuration::Mesh>(topo, conf, sim, r, size);

  DEBUG(10,"size domain(" << rank << "): " << conf.lattice_sum().domain.size());
  
  // always calculate at the beginning or read it from file.
  if (sim.steps() == 0 && !sim.param().nonbonded.influence_function_read) {
    DEBUG(10,"\t calculating influence function");
    if (rank == 0)
      timer.start("P3M: influence function");
    if (sim.mpi)
      conf.lattice_sum().influence_function.calculate<configuration::ParallelMesh>(topo, conf, sim);
    else
      conf.lattice_sum().influence_function.calculate<configuration::Mesh>(topo, conf, sim);

    const double new_rms_force_error = sqrt(conf.lattice_sum().influence_function.quality()
            / (math::volume(conf.current().box, conf.boundary_type) * topo.num_atoms()));
    if (rank == 0)
      timer.stop("P3M: influence function");
    
    if (new_rms_force_error > sim.param().nonbonded.influence_function_rms_force_error) {
      io::messages.add("P3M: RMS force error is still too big after reevaluation "
              "of the influence function. Increase the number of grid points, "
              "the rms force error threshold, or the charge width parameter.",
              "P3M", io::message::error);
      return;
    }
  }

  // check whether we have to update the influence function
  if (sim.steps() && sim.param().nonbonded.accuracy_evaluation &&
      sim.steps() % sim.param().nonbonded.accuracy_evaluation == 0) {
    if (rank == 0)
      timer.start("P3M: accuracy evaluation");
    if (sim.mpi)
      conf.lattice_sum().influence_function.evaluate_quality<configuration::ParallelMesh>(topo, conf, sim);
    else
      conf.lattice_sum().influence_function.evaluate_quality<configuration::Mesh>(topo, conf, sim);
    // number of charges is set to number of atoms as in promd...
    // see MD02.10 eq. C7
    const double rms_force_error = sqrt(conf.lattice_sum().influence_function.quality()
            / (math::volume(conf.current().box, conf.boundary_type) * topo.num_atoms()));
    if (rank == 0)
      timer.stop("P3M: accuracy evaluation");
    
    if (rms_force_error > sim.param().nonbonded.influence_function_rms_force_error) {
      // recalculate the influence function
      DEBUG(10,"\t calculating influence function");
      if (rank == 0)
        timer.start("P3M: influence function");
      if (sim.mpi)
        conf.lattice_sum().influence_function.calculate<configuration::ParallelMesh>(topo, conf, sim);
      else
        conf.lattice_sum().influence_function.calculate<configuration::Mesh>(topo, conf, sim);
      
      const double new_rms_force_error = sqrt(conf.lattice_sum().influence_function.quality()
            / (math::volume(conf.current().box, conf.boundary_type) * topo.num_atoms()));
      if (rank == 0)
        timer.stop("P3M: influence function");
      
      if (new_rms_force_error > sim.param().nonbonded.influence_function_rms_force_error) {
        io::messages.add("P3M: RMS force error is still too big after reevaluation "
                "of the influence function. Increase the number of grid points, "
                "the rms force error threshold, or the charge width parameter.",
                "P3M", io::message::error);
        return;
      }
    }
  }
  if (rank == 0)
    timer.start("P3M: energy & force");
  
  DEBUG(10,"\t done with influence function, starting to assign charge density to grid ... ");
  if (sim.mpi)
    interaction::Lattice_Sum::calculate_charge_density<configuration::ParallelMesh>(topo, conf, sim, r);
  else
    interaction::Lattice_Sum::calculate_charge_density<configuration::Mesh>(topo, conf, sim, r);

  DEBUG(10,"\t assigned charge density to grid, starting fft of charge density");
  // FFT the charge density grid
  configuration::Mesh & charge_density = *conf.lattice_sum().charge_density;
  charge_density.fft(configuration::Mesh::fft_forward);
  DEBUG(10, "\t done with fft! Starting to calculate the energy ...");
  
  if (sim.mpi) {
    interaction::Lattice_Sum::calculate_potential_and_energy<configuration::ParallelMesh>(topo, conf, sim, storage);
  } else {
    interaction::Lattice_Sum::calculate_potential_and_energy<configuration::Mesh>(topo, conf, sim, storage);
  }
  DEBUG(10, "\t done with calculation of elec. potential and energy, calculating electric field...");
  if (sim.mpi) {
    interaction::Lattice_Sum::calculate_electric_field<configuration::ParallelMesh>(topo, conf, sim);
  } else {
    interaction::Lattice_Sum::calculate_electric_field<configuration::Mesh>(topo, conf, sim);
  }
  DEBUG(10, "\t done with electric field calculation, calculating forces");
  if (sim.mpi) {
    interaction::Lattice_Sum::calculate_force<configuration::ParallelMesh>(topo, conf, sim, storage, r);
  } else {
    interaction::Lattice_Sum::calculate_force<configuration::Mesh>(topo, conf, sim, storage, r);
  }
  
  if (rank == 0)
    timer.stop("P3M: energy & force");
  DEBUG(7, "\tdone with calculating interactions in k-space (P3M)");  
}

void interaction::Nonbonded_Outerloop
::ls_self_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size) {
  DEBUG(8, "\telectrostatic self energy");
  double a1 = 0.0, a3 = 0.0;
  double & a2_tilde = conf.lattice_sum().a2_tilde;
  double a2;

  const int shape = sim.param().nonbonded.ls_charge_shape;
  // see MD05.32 Table 6
  switch (shape) {
    case -1 :
      a1 = 1.0;
      a3 = 2.0 * 1.0 / sqrt(math::Pi);
      break;
    case 0:
      a1 = 2.0 / 5.0;
      a3 = 3.0 / 2.0;
      break;
    case 1:
      a1 = 4.0 / 15.0;
      a3 = 2.0;
      break;
    case 2:
      a1 = 4.0 / 21.0;
      a3 = 5.0 / 2.0;
      break;
    case 3:
      a1 = 3.0 / 14.0;
      a3 = 9.0 / 4.0;
      break;
    case 4:
      a1 = 1.0 / 6.0;
      a3 = 21.0 / 8.0;
      break;
    case 5:
      a1 = 2.0 / 15.0;
      a3 = 3.0;
      break;
    case 6:
      a1 = 8.0 / 55.0;
      a3 = 45.0 / 16.0;
      break;
    case 7:
      a1 = 4.0 / 33.0;
      a3 = 25.0 / 8.0;
      break;
    case 8:
      a1 = 4.0 / 39.0;
      a3 = 55.0 / 16.0;
      break;
    case 9:
      a1 = 10.0 / 91.0;
      a3 = 105.0 / 32.0;
      break;
    case 10:
      a1 = 2.0 / 21.0;
      a3 = 455.0 / 128.0;
      break;
    default:
      io::messages.add("charge shape not implemented", "Lattice Sum",
              io::message::critical);
  }
  const double volume = math::volume(conf.current().box, conf.boundary_type);
  const double width = sim.param().nonbonded.ls_charge_shape_width; 
  // see MD05.32 Table 6
  a1 *= -math::Pi * width * width / volume;
  a3 *= -1.0 / width;

  // calculate the a2 term
  // do we have to do it numerically?
  if (sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2_numerical ||
      sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact_a2_numerical ||
      sim.param().nonbonded.ls_calculate_a2 == simulation::la_a2t_ave_a2_numerical) {
    // calculate A2 numerically
    const double & required_precision = sim.param().nonbonded.ls_a2_tolerance;
    math::Matrix l_to_k = configuration::KSpace_Utils::l_to_k_matrix(
            conf.current().box, conf.boundary_type);

    // Here, we loop over the surface of triclinic volumes of increasing
    // size in l-space, and add the successive A2 contributions of these
    // surfaces.
    math::GenericVec<int> l(0);
    math::Vec k;
    int l_max = 0;
    double tolerance;
    a2 = 0.0;
    do {
      ++l_max;
      double term = 0.0;

      // the planes are located perpendicular to the axis (coord)
      for (unsigned int coord = 0; coord < 3; ++coord) {
        const unsigned int coord_a = (coord + 1) % 3;
        const unsigned int coord_b = (coord + 2) % 3;
        const int boundary_a = (coord > 0) ? l_max - 1 : l_max;
        const int boundary_b = (coord > 1) ? l_max - 1 : l_max;

        // the plane can be located at -l_max or +l_max
        for (int sign = -1; sign <= 1; sign += 2) {
          DEBUG(12, "\tnew plane");
          l(coord) = sign * l_max;

          // loop over the plane excluding edges for some axes
          for (int l_a = -boundary_a; l_a <= boundary_a; ++l_a) {
            l(coord_a) = l_a;
            for (int l_b = -boundary_b; l_b <= boundary_b; ++l_b) {
              l(coord_b) = l_b;
              DEBUG(13, "\t\tl: " << math::v2s(l));
              k = math::product(l_to_k, l);
              double k2 = math::abs2(k);
              double gamma_hat;
              interaction::Lattice_Sum::charge_shape_fourier(shape,
                      sqrt(k2) * width, gamma_hat);
              term += gamma_hat / k2;
              DEBUG(13, "\t\t\tgamma_hat / k2: " << gamma_hat / k2);
            }
          } // loop over planes         
        } // loop over signs
      } // loop over coordinates

      // now calculate A2 and the relative tolerance
      a2 += term;
      tolerance = fabs(term / a2);
      DEBUG(11, "\ttolerance: " << tolerance);
    } while (tolerance > required_precision);

    a2 *= 4.0 * math::Pi / volume;
  }
  
  switch (sim.param().nonbonded.ls_calculate_a2) {
    case simulation::ls_a2_zero :
      // we just set both A2 to zero.
      a2 = a2_tilde = 0.0;
      break;
    case simulation::ls_a2t_exact :
      // A2t was calculated exactly by Ewald/P3M. A2 is set to A2t
      a2 = a2_tilde;
    case simulation::ls_a2_numerical :
      // A2 was calculate numerically, A2t is set to A2
      a2_tilde = a2;
      break;
    case simulation::ls_a2t_exact_a2_numerical :
    case simulation::la_a2t_ave_a2_numerical :
      // we already have A2t and A2 - do nothing
      break;
    default :
      io::messages.add("A2 calculation method not implemented", "Lattice Sum",
              io::message::critical);
  } // switch ls_calculate_a2
  
  // calculate box overall square charge (MD05.32 eq 41)
  const unsigned int num_atoms = topo.num_atoms();
  double st2 = 0.0, s = 0;
  for(unsigned int i = 0; i < num_atoms; ++i) {
    const double qi = topo.charge(i);
    s += qi;
    st2 += qi * qi;
  }
  DEBUG(10,"a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3);
  // now combine everything (MD05.32 eq 54)
  storage.energies.ls_self_total = (a1 + a2 + a3) * st2 * math::eps0_i / (8.0 * math::Pi);
  DEBUG(8,"ls_self_total = " << storage.energies.ls_self_total);
  
  // now claculate the E_A term
  storage.energies.ls_a_term_total = (a1 * s * s - (a1 + a2_tilde) * st2) * math::eps0_i / (8.0 * math::Pi);
  DEBUG(8, "ls_a_term_total = " << storage.energies.ls_a_term_total);
}

void interaction::Nonbonded_Outerloop
::ls_surface_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size) {
  SPLIT_INNERLOOP(_ls_surface_outerloop, topo, conf, sim,
          storage, rank, size);
}
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_ls_surface_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size)
{  
  DEBUG(8, "\telectrostatic surface terms");

  // under tinfoil boundary conditions we don't have
  // to calculate anything as the terms are zero.
  if (sim.param().nonbonded.ls_epsilon == 0.0) {
     DEBUG(10,"\ttinfoil.");
     storage.energies.ls_surface_total = 0.0;
     return;
  }

  math::Vec box_dipole_moment(0.0);
  math::Vec box_centre = conf.current().box(0) / 2.0 + 
          conf.current().box(1) / 2.0 +
          conf.current().box(2) / 2.0;

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  
  const unsigned int num_atoms = topo.num_atoms();
  for(unsigned int i = 0; i < num_atoms; i++) {
    math::Vec r = conf.current().pos(i);
    //periodicity.put_into_positive_box(r);

    box_dipole_moment += topo.charge(i) * (r - box_centre);
  }

  // now we have the box dipole moment and can calculate the energy and force
  const double prefactor = math::eps0_i /
           ((sim.param().nonbonded.ls_epsilon * 2.0 + 1.0) *
            math::volume(conf.current().box, conf.boundary_type));
 
  storage.energies.ls_surface_total = 0.5 * math::abs2(box_dipole_moment)*prefactor;
  DEBUG(10, "\tsurface energy: " << storage.energies.ls_surface_total);

  for(unsigned int i = 0; i < num_atoms; i++) {
    storage.force(i) += - prefactor * topo.charge(i) * box_dipole_moment;
  }
}

/**
 * calculate the interaction for a given atom pair.
 * SLOW! as it has to create the periodicity...
 */
int interaction::Nonbonded_Outerloop::calculate_interaction
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 unsigned int atom_i, unsigned int atom_j,
 math::Vec & force, 
 double &e_lj, double &e_crf
 )
{
  SPLIT_INNERLOOP(_calculate_interaction, topo, conf, sim, atom_i, atom_j, force, e_lj, e_crf);
  return 0;
}

template<typename t_interaction_spec>
int interaction::Nonbonded_Outerloop
::_calculate_interaction(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 unsigned int atom_i, unsigned int atom_j,
			 math::Vec & force,
			 double & e_lj, double & e_crf)
{
  math::Vec r;
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  
  Nonbonded_Term term;
  term.init(sim);

  const lj_parameter_struct &lj = 
    m_param.lj_parameter(topo.iac(atom_i),
			 topo.iac(atom_j));
  
  periodicity.nearest_image(conf.current().pos(atom_i), conf.current().pos(atom_j), r);

  double f;
  term.lj_crf_interaction(r, lj.c6, lj.c12,
			  topo.charge()(atom_i) * topo.charge()(atom_j),
			  f, e_lj, e_crf);
  force = f * r;

  return 0;
}


/**
 * calculate the hessian for a given atom.
 * this will be VERY SLOW !
 */
int interaction::Nonbonded_Outerloop
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    unsigned int atom_i, unsigned int atom_j,
		    math::Matrix & hessian,
		    PairlistContainer const & pairlist){

  SPLIT_INNERLOOP(_calculate_hessian, topo, conf, sim, atom_i, atom_j, hessian, pairlist);
  return 0;
}

template<typename t_interaction_spec>
int interaction::Nonbonded_Outerloop
::_calculate_hessian(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     unsigned int atom_i, unsigned int atom_j,
		     math::Matrix & hessian,
		     PairlistContainer const & pairlist)
{
  
  hessian = 0.0;
  
  // loop over the pairlist

  //*************************
  // standard implementation
  //*************************
  
  // check whether the pair is in one of the shortrange pairlists
  bool calculate_pair = 
    std::find(pairlist.solute_short[atom_i].begin(),
                pairlist.solute_short[atom_i].end(),
                atom_j) != pairlist.solute_short[atom_i].end() || // i-j in solute
    std::find(pairlist.solute_short[atom_j].begin(),
                pairlist.solute_short[atom_j].end(),
                atom_i) != pairlist.solute_short[atom_j].end() || // j-i in solute
    std::find(pairlist.solvent_short[atom_i].begin(),
                pairlist.solvent_short[atom_i].end(),
                atom_j) != pairlist.solvent_short[atom_i].end() || // i-j in solvent
    std::find(pairlist.solvent_short[atom_j].begin(),
                pairlist.solvent_short[atom_j].end(),
                atom_i) != pairlist.solvent_short[atom_j].end(); // j-i in solvent

  if (calculate_pair) {
    math::Vec r;
    math::Matrix h;
    math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
    
    Nonbonded_Term term;
    term.init(sim);
    
    DEBUG(8, "\thessian pair in pairlist: " << atom_i << " - " << atom_j);
    
    periodicity.nearest_image(conf.current().pos(atom_i),
            conf.current().pos(atom_j),
            r);
    
    const lj_parameter_struct &lj =
            m_param.lj_parameter(topo.iac(atom_i),
            topo.iac(atom_j));
    
    term.lj_crf_hessian(r,
            lj.c6, lj.c12,
            topo.charge()(atom_i) * topo.charge()(atom_j),
            h);
    
    for(unsigned int d1=0; d1 < 3; ++d1)
      for(unsigned int d2=0; d2 < 3; ++d2)
        hessian(d1, d2) += h(d1, d2);
  }
  
  return 0;
}
