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

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
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
    ++cg_it;

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
