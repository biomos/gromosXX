/**
 * @file perturbed_shake.cc
 * contains the template methods for
 * the class Perturbed_Shake.
 */

#ifdef XXMPI
#include <mpi.h>
#endif

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/shake.h"
#include "../../algorithm/constraints/perturbed_shake.h"

#include "../../util/template_split.h"
#include "../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

#include "perturbed_dihedral_constraint.cc"

/**
 * Constructor.
 */
algorithm::Perturbed_Shake::Perturbed_Shake(double const tolerance, int const max_iterations)
  : Shake(tolerance, max_iterations)
{
}

/**
 * Destructor.
 */
algorithm::Perturbed_Shake::~Perturbed_Shake()
{
}

//================================================================================
//         PERTURBED SHAKE ITERATION
//================================================================================

/**
 * do one iteration
 */      
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Perturbed_Shake
::perturbed_shake_iteration(topology::Topology const &topo,
			    configuration::Configuration & conf,
			    bool & convergence,
			    int const first,
			    std::vector<bool> &skip_now,
			    std::vector<bool> &skip_next,
			    std::vector<topology::perturbed_two_body_term_struct>
			    const & constr,
			    double const dt,
			    math::Periodicity<B> const & periodicity,
			    simulation::Simulation & sim) // ANITA
{
  convergence = true;

  // index for constraint_force...
  unsigned int k = 0;
  double const dt2 = dt * dt;
  
  // and constraints
  for(typename std::vector<topology::perturbed_two_body_term_struct>
	::const_iterator
	it = constr.begin(),
	to = constr.end();
      it != to;
      ++it, ++k ){
	
    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j]) continue;
    if (topo.inverse_mass(first + it->i) == 0 &&
            topo.inverse_mass(first + it->j) == 0) continue;

    DEBUG(10, "i: " << it->i << " j: " << it->j << " first: " << first);

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    // we use the bond lambda
    const double lam = topo.individual_lambda(simulation::bond_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];
    const double lam_derivative = topo.individual_lambda_derivative
      (simulation::bond_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];

    // the position
    math::Vec &pos_i = conf.current().pos(first+it->i);
    math::Vec &pos_j = conf.current().pos(first+it->j);

    DEBUG(10, "\ni: " << math::v2s(pos_i) << "\nj: " << math::v2s(pos_j));
	
    math::Vec r;
    periodicity.nearest_image(pos_i, pos_j, r);
    double dist2 = math::abs2(r);
	
    double r0 = (1.0 - lam) * topo.bond_types_harm()[it->A_type].r0 + 
      lam * topo.bond_types_harm()[it->B_type].r0;

    DEBUG(10, "constraint length: " << r0);
    DEBUG(10, "r0(A) = " << topo.bond_types_harm()[it->A_type].r0);
    DEBUG(10, "r0(B) = " << topo.bond_types_harm()[it->B_type].r0);    

    double constr_length2 = r0 * r0;
    double diff = constr_length2 - dist2;

    DEBUG(15, "constr: " << constr_length2 << " dist2: " << dist2);
	  
    if(fabs(diff) >= constr_length2 * this->tolerance() * 2.0){
      // we have to shake
      DEBUG(10, "shaking");
      
      // the reference position
      const unsigned int atom_i = first+it->i;
      const unsigned int atom_j = first+it->j;
      const math::Vec &ref_i = conf.old().pos(atom_i);
      const math::Vec &ref_j = conf.old().pos(atom_j);
    DEBUG(10, "\niref: " << math::v2s(ref_i) << "\njref: " << math::v2s(ref_j));
      
      math::Vec ref_r;
      periodicity.nearest_image(ref_i, ref_j, ref_r);

      double sp = dot(ref_r, r);
	  
      if(sp < constr_length2 * math::epsilon){
	io::messages.add("SHAKE error. vectors orthogonal",
			 "Shake::???",
			 io::message::critical);
	DEBUG(5, "ref i " << math::v2s(ref_i) << " ref j " << math::v2s(ref_j));
	DEBUG(5, "free i " << math::v2s(pos_i) << " free j " << math::v2s(pos_j));
	DEBUG(5, "ref r " << math::v2s(ref_r));
	DEBUG(5, "r " << math::v2s(r));
	
	std::cout << "Perturbed SHAKE ERROR\n"
		  << "\tatom i    : " << atom_i + 1 << "\n"
		  << "\tatom j    : " << atom_j + 1 << "\n"
		  << "\tfirst     : " << first << "\n"
		  << "\tref i     : " << math::v2s(ref_i) << "\n"
		  << "\tref j     : " << math::v2s(ref_j) << "\n"
		  << "\tfree i    : " << math::v2s(pos_i) << "\n"
		  << "\tfree j    : " << math::v2s(pos_j) << "\n"
		  << "\tref r     : " << math::v2s(ref_r) << "\n"
		  << "\tr         : " << math::v2s(r) << "\n"
		  << "\tsp        : " << sp << "\n"
		  << "\tconstr    : " << constr_length2 << "\n"
		  << "\tdiff      : " << diff << "\n"
		  << "\tforce i   : " << math::v2s(conf.old().force(atom_i)) << "\n"
		  << "\tforce j   : " << math::v2s(conf.old().force(atom_j)) << "\n"
		  << "\tvel i     : " << math::v2s(conf.current().vel(atom_i)) << "\n"
		  << "\tvel j     : " << math::v2s(conf.current().vel(atom_j)) << "\n"
		  << "\told vel i : " << math::v2s(conf.old().vel(atom_i)) << "\n"
		  << "\told vel j : " << math::v2s(conf.old().vel(atom_j)) << "\n\n";
	
	return E_SHAKE_FAILURE;
      }
	  
      // lagrange multiplier
      double lambda = diff / (sp * 2 *
			      (topo.inverse_mass()(atom_i) +
			       topo.inverse_mass()(atom_j) ));      

      DEBUG(10, "lagrange multiplier " << lambda);

      const math::Vec cons_force = lambda * ref_r;
      conf.old().constraint_force(atom_i) += cons_force;
      conf.old().constraint_force(atom_j) -= cons_force;

      if (V == math::atomic_virial) {
        for (int a = 0; a < 3; ++a) {
          for (int aa = 0; aa < 3; ++aa) {
            conf.old().virial_tensor(a, aa) +=
                    ref_r(a) * ref_r(aa) * lambda / dt2;
          }
        }
        DEBUG(12, "\tatomic virial done");
      }
      
      // the perturbed energy derivatives

      conf.old().perturbed_energy_derivatives.
	constraints_energy[topo.atom_energy_group()[it->i]] +=
	lam_derivative * lambda / dt2 * sqrt(constr_length2) *
	(topo.bond_types_harm()[it->B_type].r0 - topo.bond_types_harm()[it->A_type].r0);

      //ANITA
      if (sim.param().precalclam.nr_lambdas &&
          ((sim.steps()  % sim.param().write.free_energy) == 0)){
        DEBUG(1, "AB_bond within if ");
        double r0A = topo.bond_types_harm()[it->A_type].r0;
        DEBUG(1, "AB_bond r0A " << r0A); 
        double r0B = topo.bond_types_harm()[it->B_type].r0;
        DEBUG(1, "AB_bond r0B " << r0B); 
        double r0_diff = r0B - r0A; 
        DEBUG(1, "AB_bond r0_diff " << r0_diff); 
        double sp_2_m_dt2 = sp * 2 * (1.0 / topo.mass()(atom_i) +
                               1.0 / topo.mass()(atom_j) ) * dt2; 
        DEBUG(1, "AB_bond sp_2_m_dt2 " << sp_2_m_dt2); 

        double lambda_step = (sim.param().precalclam.max_lam -
                   sim.param().precalclam.min_lam) /
                   (sim.param().precalclam.nr_lambdas-1);

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

          // determine current lambda for this index
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

          double r0lam = (1-lam)*r0A + lam*r0B;
          DEBUG(1, "AB_bond r0lam " << r0lam); 
          double difflam = r0lam*r0lam - dist2; 
          DEBUG(1, "AB_bond difflam " << difflam); 
          double value = (difflam / sp_2_m_dt2) * r0lam *r0_diff;
          conf.old().perturbed_energy_derivatives.AB_bond[lam_index] += value;
        }
      } // ANITA 

      // update positions
      ref_r *= lambda;
      pos_i += ref_r * topo.inverse_mass()(first+it->i);
      pos_j -= ref_r * topo.inverse_mass()(first+it->j);
	  
      convergence = false;

      // consider atoms in the next step
      skip_next[it->i] = false;
      skip_next[it->j] = false;
      
    } // we have to shake
  } // constraints
  
  return 0;
}    

//================================================================================
//         PERTURBED SHAKE SOLUTE / SOLVENT LOOPS
//================================================================================

/**
 * shake perturbed solute
 */
template<math::boundary_enum B, math::virial_enum V>
void algorithm::Perturbed_Shake
::perturbed_solute(topology::Topology const & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
//ANITA		   simulation::Simulation const & sim,
		   int max_iterations,
		   int &error)
{
  // for now shake the whole solute in one go,
  // not bothering about submolecules...

  if (!sim.mpi || m_rank == 0)
    m_timer.start("solute");

  DEBUG(8, "\tshaking perturbed SOLUTE"); 
  math::Periodicity<B> periodicity(conf.current().box);
  
  const unsigned int num_atoms = topo.num_solute_atoms();   
  std::vector<bool> skip_now;
  std::vector<bool> skip_next;
  skip_next.assign(topo.solute().num_atoms(), true);
  skip_now.assign(topo.solute().num_atoms(), false);

  int first = 0;
  error = 0;
  
  int num_iterations = 0;
  bool convergence = false;
  while(!convergence){
    DEBUG(9, "\titeration" << std::setw(10) << num_iterations);

    // distance constraints
    bool dist_convergence = true, pert_dist_convergence = true;
    
    if (topo.perturbed_solute().distance_constraints().size() && 
	sim.param().constraint.solute.algorithm == simulation::constr_shake &&
	sim.param().constraint.ntc > 1){

      DEBUG(8, "perturbed shake iteration (solute distance)");
      if(perturbed_shake_iteration<B, V>
	 (topo, conf, pert_dist_convergence, first, skip_now, skip_next,
	  topo.perturbed_solute().distance_constraints(), sim.time_step_size(),
	  periodicity, sim)){
	io::messages.add("Perturbed SHAKE error. vectors orthogonal",
			 "Perturbed_Shake::solute",
			 io::message::error);
	std::cout << "Perturbed SHAKE failure in solute!" << std::endl;
	error = E_SHAKE_FAILURE_SOLUTE;
	break;
      }
    }
    
    if (topo.solute().distance_constraints().size() && 
	sim.param().constraint.solute.algorithm == simulation::constr_shake &&
	sim.param().constraint.ntc > 1){
      
      DEBUG(8, "unperturbed shake iteration (solute distance)");
      if(Shake::shake_iteration<B, V>
	 (topo, conf, dist_convergence, 
	  first, skip_now, skip_next,
	  topo.solute().distance_constraints(), sim.time_step_size(),
	  periodicity) != 0){

	io::messages.add("SHAKE error. vectors orthogonal",
			 "Shake::solute",
			 io::message::error);
	std::cout << "SHAKE failure in solute!" << std::endl;
	error = E_SHAKE_FAILURE_SOLUTE;
	break;
      }
    }
    
    // dihedral constraints
    bool dih_convergence = true, pert_dih_convergence = true;
    if (sim.param().dihrest.dihrest == simulation::dihedral_constr) {      
      DEBUG(7, "SHAKE: perturbed dihedral constraints iteration");
      if(perturbed_dih_constr_iteration<B, V>
	 (topo, conf, sim, pert_dih_convergence, skip_now, skip_next, periodicity)
	 ){
	io::messages.add("SHAKE error: perturbed dihedral constraints",
			 "Shake::solute",
			 io::message::error);
	std::cout << "SHAKE failure in solute perturbed dihedral constraints!" << std::endl;
	error = E_SHAKE_FAILURE_SOLUTE;
	break;
      }

      DEBUG(7, "SHAKE: dihedral constraints iteration");
      if(dih_constr_iteration<B, V>
	 (topo, conf, sim, dih_convergence, skip_now, skip_next,topo.dihedral_restraints(), periodicity)
	 ){
	io::messages.add("SHAKE error: dihedral constraints",
			 "Shake::solute",
			 io::message::error);
	std::cout << "SHAKE failure in solute dihedral constraints!" << std::endl;
	error = E_SHAKE_FAILURE_SOLUTE;
	break;
      }
    }

    convergence = pert_dist_convergence && dist_convergence 
      && pert_dih_convergence && dih_convergence;

    if(++num_iterations > max_iterations){
      io::messages.add("Perturbed SHAKE error. too many iterations",
		       "Perturbed_Shake::solute",
		       io::message::error);
      error = E_SHAKE_FAILURE_SOLUTE;
      break;
    }

    skip_now = skip_next;
    skip_next.assign(skip_next.size(), true);

  } // convergence?
  
  // constraint force
  const double dt2 = sim.time_step_size() * sim.time_step_size();
  for (unsigned int i=0; i < num_atoms; ++i){
    conf.old().constraint_force(i) *= 1 / dt2;
    DEBUG(5, "constraint_force " << math::v2s(conf.old().constraint_force(i)));
  }

  if (!sim.mpi || m_rank == 0)
    m_timer.stop("solute");

  //error = 0;

} // solute

//================================================================================
//         apply PERTURBED SHAKE
//================================================================================

/**
 * apply the SHAKE algorithm
 */
int algorithm::Perturbed_Shake
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying perturbed SHAKE");
  
  if (!sim.mpi || m_rank == 0)
    m_timer.start();
  
  // set the constraint force to zero
  std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
  for (; it != to; ++it) {
    conf.old().constraint_force(*it) = 0.0;
  }
  
  int error = 0;
  
  if (m_rank == 0) {
    // check whether we shake solute
    if (((topo.perturbed_solute().distance_constraints().size() ||
            topo.solute().distance_constraints().size()) &&
            sim.param().constraint.solute.algorithm == simulation::constr_shake &&
            sim.param().constraint.ntc > 1) ||
            sim.param().dihrest.dihrest == simulation::dihedral_constr) {
      
      DEBUG(8, "\twe need to shake perturbed SOLUTE");
      
      SPLIT_VIRIAL_BOUNDARY(perturbed_solute,
              topo, conf, sim,
              this->max_iterations(), error);
    }
  }
  
  // broadcast eventual errors from master to slaves
#ifdef XXMPI
  if (sim.mpi) {
    MPI_Bcast(&error, 1, MPI::INT, sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
  } 
#endif

  if (error){
      std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLUTE "
      << "at step " << sim.steps() << std::endl;
      conf.special().shake_failure_occurred = true;
      m_timer.stop();
      return E_SHAKE_FAILURE_SOLUTE;
  }
  

#ifdef XXMPI
  math::VArray & pos = conf.current().pos;
  if (sim.mpi) {
    // broadcast current and old coordinates and pos.

    MPI_Bcast(&pos(0)(0), pos.size() * 3, MPI::DOUBLE, sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
    MPI_Bcast(&conf.old().pos(0)(0), conf.old().pos.size() * 3, MPI::DOUBLE, sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
    MPI_Bcast(&conf.current().box(0)(0), 9, MPI::DOUBLE, sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);

    // set virial tensor and solute coordinates of slaves to zero
    if (m_rank) { // slave
      conf.old().virial_tensor = 0.0;
      conf.old().constraint_force = 0.0;
      
    for(unsigned int i = 0; i <topo.num_solute_atoms(); ++i){
      pos(i) = 0;
    }
    }
  }
#endif

  // solvent
  bool do_shake_solv=sim.param().system.nsm &&
      sim.param().constraint.solvent.algorithm == simulation::constr_shake;
  if (do_shake_solv){

    DEBUG(8, "\twe need to shake SOLVENT");
    
    SPLIT_VIRIAL_BOUNDARY(solvent, 
			  topo, conf, sim, sim.time_step_size(), 
			  this->max_iterations(), error);
  }
  
  if (error){
    std::cout << "SHAKE: exiting with error condition: E_SHAKE_FAILURE_SOLVENT "
	      << "at step " << sim.steps() << std::endl;
    conf.special().shake_failure_occurred = true;
    return E_SHAKE_FAILURE_SOLVENT;
  }
#ifdef XXMPI
  if (sim.mpi && do_shake_solv) {
    if (m_rank == 0) {
      // Master 
      // reduce current positions, store them in new_pos and assign them to current positions
      math::VArray new_pos=pos;
      MPI_Reduce(&pos(0)(0), &new_pos(0)(0),
              topo.num_atoms() * 3, MPI::DOUBLE, MPI::SUM,  sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
      pos = new_pos;

      // reduce current virial tensor, store it in virial_new and reduce it to current tensor
      math::Matrix virial_new(0.0);
      MPI_Reduce(&conf.old().virial_tensor(0, 0), &virial_new(0, 0),
              9, MPI::DOUBLE, MPI::SUM,  sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
      conf.old().virial_tensor = virial_new;

      // reduce current contraint force, store it in cons_force_new and reduce
      // it to the current constraint force
      math::VArray cons_force_new=conf.old().constraint_force;
      MPI_Reduce(&conf.old().constraint_force(0)(0), &cons_force_new(0)(0),
              topo.num_atoms() * 3, MPI::DOUBLE, MPI::SUM,  sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
      conf.old().constraint_force = cons_force_new;
    } else {
      // slave
      // reduce pos
      MPI_Reduce(&pos(0)(0), NULL,
              topo.num_atoms() * 3, MPI::DOUBLE, MPI::SUM, sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
      // reduce virial
      MPI_Reduce(&conf.old().virial_tensor(0, 0), NULL,
              9, MPI::DOUBLE, MPI::SUM, sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
      // reduce constraint force
      MPI_Reduce(&conf.old().constraint_force(0)(0), NULL,
              topo.num_atoms() * 3, MPI::DOUBLE, MPI::SUM, sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);
    }
  }
#endif
  
  // shaken velocity:
  // stochastic dynamics, energy minimisation, analysis needs to shake without
  // velocity correction (once; it shakes twice...)
  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze) {
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.current().vel(*it) = (conf.current().pos(*it) - conf.old().pos(*it)) /
              sim.time_step_size();
    }
  }
  
  // return success!
  if (!sim.mpi || m_rank == 0)
    m_timer.stop();
  return error;	   
}

//================================================================================
//         PERTURBED SHAKE INITIALIZATION
//================================================================================

int algorithm::Perturbed_Shake::init(topology::Topology & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation & sim,
				     std::ostream & os,
				     bool quiet)
{
  if (!quiet){
    os << "PERTURBED SHAKE\n"
	    << "\tsolute\t";
    if (sim.param().constraint.solute.algorithm == simulation::constr_shake
	&& topo.perturbed_solute().distance_constraints().size()){    
      os << "ON\n";  
      os << "\t\ttolerance = "
		<< sim.param().constraint.solute.shake_tolerance << "\n";
    }
    else os << "OFF\n";
  }
  
  #ifdef XXMPI
  if (sim.mpi) {
    m_rank = sim.mpi_control.simulationThisThreadID;
    m_size = sim.mpi_control.simulationNumberOfThreads;
  } else {
    m_rank = 0;
    m_size = 1;
  }
#else
  m_rank = 0;
  m_size = 1;
#endif

  if (sim.param().constraint.solute.algorithm == simulation::constr_shake) {
    // loop over the constraints to find out which atoms are constrained
    {
      std::vector<topology::two_body_term_struct>::const_iterator
      it = topo.solute().distance_constraints().begin(),
              to = topo.solute().distance_constraints().end();
      for (; it != to; ++it) {
        constrained_atoms().insert(it->i);
        constrained_atoms().insert(it->j);
      }
    }
    // loop over the perturbed constraints
    {
      std::vector<topology::perturbed_two_body_term_struct>::const_iterator
      it = topo.perturbed_solute().distance_constraints().begin(),
              to = topo.perturbed_solute().distance_constraints().end();
      for (; it != to; ++it) {
        constrained_atoms().insert(it->i);
        constrained_atoms().insert(it->j);
      }
    }   
    // also add the dihedral constrained atoms
    if (sim.param().dihrest.dihrest == simulation::dihedral_constr) {
      {
        std::vector<topology::dihedral_restraint_struct>::const_iterator
        it = topo.dihedral_restraints().begin(),
                to = topo.dihedral_restraints().end();
        for (; it != to; ++it) {
          constrained_atoms().insert(it->i);
          constrained_atoms().insert(it->j);
          constrained_atoms().insert(it->k);
          constrained_atoms().insert(it->l);
        }
      }
      { // and perturbed dihedrals
        std::vector<topology::perturbed_dihedral_restraint_struct>::const_iterator
        it = topo.perturbed_dihedral_restraints().begin(),
                to = topo.perturbed_dihedral_restraints().end();
        for (; it != to; ++it) {
          constrained_atoms().insert(it->i);
          constrained_atoms().insert(it->j);
          constrained_atoms().insert(it->k);
          constrained_atoms().insert(it->l);
        }
      }
    }
  }
  
  if (sim.param().constraint.solvent.algorithm == simulation::constr_shake) {
    // this means that all solvent atoms are constrained. Add the to the list
    for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
      constrained_atoms().insert(i);
    }
  }
  
  if (sim.param().start.shake_pos) {
    if (!quiet)
      os << "\n\tshaking initial positions\n";

    // old and current pos and vel are the same for constrained atoms...
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.old().pos(*it) = conf.current().pos(*it);
      conf.old().vel(*it) = conf.current().vel(*it);
    }

    // shake the current ones
    if (apply(topo, conf, sim))
      return E_SHAKE_FAILURE;

    it = constrained_atoms().begin();
    for (; it != to; ++it) {
      // restore the velocities
      conf.current().vel(*it) = conf.old().vel(*it);
      // take a step back
      conf.old().pos(*it) = conf.current().pos(*it);
    }

    if (sim.param().start.shake_vel) {
      if (!quiet)
        os << "\tshaking initial velocities\n";

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        conf.current().pos(*it) = conf.old().pos(*it) -
                sim.time_step_size() * conf.old().vel(*it);
      }

      // shake again
      if (apply(topo, conf, sim))
        return E_SHAKE_FAILURE;

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        // restore the positions
        conf.current().pos(*it) = conf.old().pos(*it);
        // velocities are in opposite direction (in time)
        conf.current().vel(*it) = -1.0 * conf.current().vel(*it);
        conf.old().vel(*it) = conf.current().vel(*it);
      }
    } // if shake vel
  } else if (sim.param().start.shake_vel) {
    io::messages.add("shaking velocities without shaking positions illegal.",
            "shake", io::message::error);
  }

  if (!quiet) {
    os << "END\n"
       << "SHAKE\n" 
       << "\tsolvent\t";
  
    if (sim.param().constraint.solvent.algorithm == simulation::constr_shake){
      if (sim.mpi)
        os << "ON (MPI parallel version)\n";
      else 
        os << "ON\n";
      os << "\t\ttolerance = " 
		<< sim.param().constraint.solvent.shake_tolerance << "\n";
    }  else os << "OFF\n";
    os << "END\n";
  }
  
  return 0;
}
