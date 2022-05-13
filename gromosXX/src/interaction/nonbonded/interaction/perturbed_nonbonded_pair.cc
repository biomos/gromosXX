/**
 * @file perturbed_nonbonded_pair.cc
 * template methods of Perturbed_Nonbonded_Pair
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/interaction_spec.h"

#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "../../../interaction/nonbonded/interaction/storage.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_innerloop.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_pair.h"

#include "../../../interaction/nonbonded/interaction_spec.h"

#include "../../../util/debug.h"
#include "../../../interaction/nonbonded/innerloop_template.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

interaction::Perturbed_Nonbonded_Pair
::Perturbed_Nonbonded_Pair(Nonbonded_Parameter &nbp)
  : m_param(&nbp)
{
}

/**
 * calculate the interactions for the
 * PERTURBED PAIRS
 * (different interaction types in A and in B)
 */
void interaction::Perturbed_Nonbonded_Pair
::perturbed_pair_outerloop(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   Storage & storage)
{
  SPLIT_PERTURBATION(_perturbed_pair_outerloop,
		     topo, conf, sim, storage);
}

template<typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Pair
::_perturbed_pair_outerloop(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_split_perturbed_pair_outerloop,
		       topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Pair
::_split_perturbed_pair_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage)
{
  DEBUG(8, "perturbed pairs");
  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  m_nonbonded_term.init(sim);
  m_perturbed_nonbonded_term.init(sim);
  // Chris:
  // This should be done when the interactions are calculated (in the innerloop, or further down if
  // we know which energy groups our atoms belong to
  // m_perturbed_nonbonded_term.set_lambda(topo.lambda(), topo.lambda_exp());

  std::vector<topology::perturbed_two_body_term_struct>::const_iterator
    it = topo.perturbed_solute().atompairs().begin(),
    to = topo.perturbed_solute().atompairs().end();
    
  for(; it != to; ++it){
    perturbed_pair_interaction_innerloop<t_interaction_spec, t_perturbation_details>
      (topo, conf, sim, it, periodicity);
  }
}

template<typename t_interaction_spec, typename perturbation_details, math::boundary_enum t_boundary_type>
void interaction::Perturbed_Nonbonded_Pair
::perturbed_pair_interaction_innerloop
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  std::vector<topology::perturbed_two_body_term_struct>
  ::const_iterator const &it,
  math::Periodicity<t_boundary_type> const & periodicity)
{

  // NO RANGE FILTER FOR PERTURBED PAIRS ??
  // NO SCALING for PERTURBED PAIRS ??
  // NO MOLECULAR VIRIAL CONTRIBUTION ??
 
  math::Vec r, f, A_f, B_f, rp1, rp2, rpp;
  math::Vec rij, rik, rjj, rjk, rim, rjm, rm;
  double A_e_lj = 0.0, A_e_crf = 0.0, A_de_lj = 0.0, A_de_crf = 0.0, 
    B_e_lj = 0.0, B_e_crf = 0.0, B_de_lj = 0.0, B_de_crf = 0.0;
  double e_lj = 0.0, e_crf = 0.0, de_lj = 0.0, de_crf = 0.0;
  lj_parameter_struct const *A_lj = nullptr;
  lj_parameter_struct const *B_lj = nullptr;
  double A_q = 0.0, B_q = 0.0, A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0;
  double alpha_lj=0, alpha_crf=0;
  
  //double A_f_pol[4], B_f_pol[4];
  math::VArray A_f_pol_vec(4), B_f_pol_vec(4), f_pol_vec(4);
  f_pol_vec = A_f_pol_vec = B_f_pol_vec = 0.0;

  
  bool is_perturbed = 0;


  DEBUG(7, "\tperturbed-pair\t" << it->i << "\t" << it->j);
  
  periodicity.nearest_image(conf.current().pos(it->i), 
			    conf.current().pos(it->j), r);

  if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func) {
    rp1 = r - conf.current().posV(it->j);
    rp2 = r + conf.current().posV(it->i);
    rpp = r + conf.current().posV(it->i) - conf.current().posV(it->j);
  }
 if (t_interaction_spec::interaction_func == simulation::pol_off_lj_crf_func) {
    periodicity.nearest_image(conf.current().pos(it->i),
            conf.current().pos(topo.gamma_j(it->i)), rij);
    periodicity.nearest_image(conf.current().pos(it->i),
            conf.current().pos(topo.gamma_k(it->i)), rik);
    periodicity.nearest_image(conf.current().pos(it->j),
            conf.current().pos(topo.gamma_j(it->j)), rjj);
    periodicity.nearest_image(conf.current().pos(it->j),
            conf.current().pos(topo.gamma_k(it->j)), rjk);

    rim = topo.gamma(it->i)*(rij + rik) / 2;
    rjm = topo.gamma(it->j)*(rjj + rjk) / 2;

    rm = r + rjm - rim;
    rp1 = rm - conf.current().posV(it->j);
    rp2 = rm + conf.current().posV(it->i);
    rpp = rm + conf.current().posV(it->i) - conf.current().posV(it->j);

  } 
  // is i perturbed?
  if (topo.is_perturbed(it->i)){
    assert(topo.perturbed_solute().atoms().count(it->i) == 1);
    
    is_perturbed = true;
      
    // is j perturbed as well?
    if(topo.is_perturbed(it->j)){
      assert(topo.perturbed_solute().atoms().count(it->j) == 1);
	
      A_lj =  &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     topo.perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].B_IAC(),
		     topo.perturbed_solute().atoms()
		     [it->j].B_IAC());

      A_qi = topo.perturbed_solute().atoms()[it->i].A_charge();
      A_qj = topo.perturbed_solute().atoms()[it->j].A_charge();

      B_qi = topo.perturbed_solute().atoms()[it->i].B_charge();
      B_qj = topo.perturbed_solute().atoms()[it->j].B_charge();
      
      A_q = A_qi * A_qj;
      B_q = B_qi * B_qj;
      
      alpha_lj = (topo.perturbed_solute().atoms()
              [it->i].LJ_softcore() +
              topo.perturbed_solute().atoms()
		  [it->j].LJ_softcore()) / 2.0;

      alpha_crf = (topo.perturbed_solute().atoms()
		   [it->i].CRF_softcore() +
		   topo.perturbed_solute().atoms()
		   [it->j].CRF_softcore() ) / 2.0;
      
    }
    else{ // only i perturbed

      A_lj = &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     topo.iac(it->j));

      B_lj = &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].B_IAC(),
		     topo.iac(it->j));

      A_qi = topo.perturbed_solute().atoms()[it->i].A_charge();
      A_qj = topo.charge()(it->j);
      
      B_qi = topo.perturbed_solute().atoms()[it->i].B_charge();
      B_qj = topo.charge()(it->j);
      
      A_q = A_qi * A_qj;
      B_q = B_qi * B_qj;
      
      alpha_lj = topo.perturbed_solute().atoms()
	[it->i].LJ_softcore();

      alpha_crf = topo.perturbed_solute().atoms()
	[it->i].CRF_softcore();
      
    }
  }
  else{ // i unperturbed
    // is j perturbed
    if(topo.is_perturbed(it->j)){
      assert(topo.perturbed_solute().atoms().count(it->j)
	       == 1);
	
      is_perturbed = true;

      A_lj =  &m_param->
	lj_parameter(topo.iac(it->i),
		     topo.perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_param->
	lj_parameter(topo.iac(it->i),
		     topo.perturbed_solute().atoms()
		     [it->j].B_IAC());

      A_qi = topo.charge()(it->i);
      A_qj = topo.perturbed_solute().atoms()[it->j].A_charge();

      B_qi = topo.charge()(it->i);
      B_qj = topo.perturbed_solute().atoms()[it->j].B_charge();
      
      A_q = A_qi * A_qj;
      B_q = B_qi * B_qj;

      alpha_lj = topo.perturbed_solute().atoms()
	[it->j].LJ_softcore();

      alpha_crf = topo.perturbed_solute().atoms()
	[it->j].CRF_softcore();
      
    }
    else{
      // both unperturbed
      
      is_perturbed = false;
	
      A_lj = &m_param->
	lj_parameter(topo.iac(it->i),
		     topo.iac(it->j));

      B_lj = A_lj;
      
      A_qi = B_qi = topo.charge()(it->i);
      A_qj = B_qj = topo.charge()(it->j);
      
      B_q = A_q = A_qi * A_qj;
    }
          
  } // all parameters done

  // calculate the lambdas
  int n1 = topo.atom_energy_group(it->i);
  int n2 = topo.atom_energy_group(it->j);
  
  // ugly check
  // we do not store the forces seperately for lj and crf for state A and B here, so we can only 
  // them with individual lambdas if the lj and crf dependencies are the same.
  if(topo.individual_lambda(simulation::lj_lambda)[n1][n2] !=
     topo.individual_lambda(simulation::crf_lambda)[n1][n2])
    
    io::messages.add("Perturbed atom pairs with different individual lambdas for lj and crf not implemented", "perturbed_nonbonded_pair", io::message::critical);
  
  // calculate all the lambda values for lj and crf interactions, softness and A and B
  m_perturbed_nonbonded_term.set_lambda
    (topo.individual_lambda(simulation::lj_lambda)[n1][n2],
     topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
     topo.individual_lambda(simulation::crf_lambda)[n1][n2],
     topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
     topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
     topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
     topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
     topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
     topo.lambda_exp());

  // interaction in state A
  // ======================
  switch(it->A_type){
    // --------------------------
    case 0: // excluded
      A_e_lj = A_e_crf = A_de_lj = A_de_crf = 0.0;
      A_f = 0.0;
      DEBUG(7, "excluded in A");
      if(sim.param().nonbonded.rf_excluded){
        if (is_perturbed){ //perturbed atoms are involved
          switch(t_interaction_spec::interaction_func){
            case simulation::lj_crf_func : {
              m_perturbed_nonbonded_term.
              rf_soft_interaction(r, A_q, 0, alpha_crf,
				  A_f, A_e_crf, A_de_crf);

      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
        double A_e_rf = 0.0, B_e_rf = 0.0, A_de_rf = 0.0, B_de_rf = 0.0;

        // determine lambda stepsize from min,max and nr of lambdas
        double lambda_step = (sim.param().precalclam.max_lam -
                 sim.param().precalclam.min_lam) /
                 (sim.param().precalclam.nr_lambdas-1);

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

          // determine current lambda for this index
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
          m_perturbed_nonbonded_term.rf_soft_interaction_ext(r, A_q, 0,
                  alpha_crf, A_e_rf, B_e_rf, 
                  A_de_rf, B_de_rf, lam);
          conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] +=  A_e_rf;
          conf.current().perturbed_energy_derivatives.A_crf_energy[lam_index]
            [topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += A_de_rf;
        }
      } // ANITA 
      
              break;
            }
            case simulation::pol_lj_crf_func : {
              /* polarization function not correctly implemented --martina
	      m_perturbed_nonbonded_term.
              pol_rf_soft_interaction(r, rp1, rp2, rpp,
                      A_qi, A_qj, B_qi, B_qj,
                      topo.coscharge(it->i), topo.coscharge(it->j),
                      alpha_crf, A_f_pol, A_e_crf, A_de_crf);
	      */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
              break;
            }
            case simulation::pol_off_lj_crf_func : {
              /* polarization function not correctly implemented --martina
	      m_perturbed_nonbonded_term.
              pol_rf_soft_interaction(rm, rp1, rp2, rpp,
                      A_qi, A_qj, B_qi, B_qj,
                      topo.coscharge(it->i), topo.coscharge(it->j),
                      alpha_crf, A_f_pol, A_e_crf, A_de_crf);
	      */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
              break;
            }
            default: io::messages.add("Perturbed_Nonbonded_Pair",
                    "interaction function not implemented",
                    io::message::critical);
          }
        } else { // non perturbed. 
          switch(t_interaction_spec::interaction_func){
            case simulation::lj_crf_func : {
              m_nonbonded_term.rf_interaction(r, A_q, A_f, A_e_crf);

	      //Since nonbonded_term is not scaled by lambda, this is calculated here --martina
	      A_f = m_perturbed_nonbonded_term.A_crf_lambda_n() * A_f;

	      A_de_crf = - topo.lambda_exp() * m_perturbed_nonbonded_term.A_crf_lambda_n_1() * A_e_crf; //define before multiplication of A_e_crf with lambda!
	      
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
         
          conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += A_e_crf;
            
        }
      } // ANITA 
      
	      A_e_crf *= m_perturbed_nonbonded_term.A_crf_lambda_n(); //scale A_e_crf with lambda

              break;
            }
            case simulation::pol_lj_crf_func : {
              /*m_nonbonded_term.
              pol_rf_interaction(r, rp1, rp2, rpp, topo.charge()(it->i),
                      topo.charge()(it->j), topo.coscharge(it->i),
                      topo.coscharge(it->j), A_f_pol_vec, A_e_crf);
	      */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
              break;
            }
            case simulation::pol_off_lj_crf_func : {
              /*m_nonbonded_term.
              pol_rf_interaction(rm, rp1, rp2, rpp, topo.charge()(it->i),
                      topo.charge()(it->j), topo.coscharge(it->i),
                      topo.coscharge(it->j), A_f_pol_vec, A_e_crf);
	      */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
              break;
            }
            default: io::messages.add("Perturbed_Nonbonded_Pair",
                    "interaction function not implemented",
                    io::message::critical);
          }
        }
      } //if(sim.param().nonbonded.rf_excluded)

      break;
      // --------------------------
    case 1: // normal interaction
      DEBUG(7, "\tlj-parameter state A c6=" << A_lj->c6 
	    << " c12=" << A_lj->c12);
      DEBUG(7, "\tcharges state A i*j = " << A_q);
      
      if (is_perturbed){
	DEBUG(7, "perturbed interaction");
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double A_f1 = 0.0, A_f6 = 0.0, A_f12 = 0.0;
            
            m_perturbed_nonbonded_term.
            lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
				    0, 0,
				    A_q, 0,
				    alpha_lj, alpha_crf,
				    A_f1, A_f6, A_f12,
				    A_e_lj, A_e_crf, A_de_lj, A_de_crf);
            
            A_f = (A_f1 + A_f6 + A_f12) * r;
            
      // ANITA
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
//        if ( sim.param().precalclam.nr_lambdas ) { 
          double A_e_lj_l = 0.0, B_e_lj_l = 0.0, A_e_crf_l = 0.0, B_e_crf_l = 0.0,
              A_de_lj_l = 0.0, B_de_lj_l = 0.0, A_de_crf_l = 0.0, B_de_crf_l = 0.0;

          // determine lambda stepsize from min,max and nr of lambdas
          double lambda_step = (sim.param().precalclam.max_lam - 
                   sim.param().precalclam.min_lam) / 
                   (sim.param().precalclam.nr_lambdas-1);

          //loop over nr_lambdas
          for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){ 

            // determine current lambda for this index
            double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

            // start the calculations
            m_perturbed_nonbonded_term.
            lj_crf_soft_interaction_ext(r, A_lj->c6, A_lj->c12,
                0, 0, A_q, 0, alpha_lj, alpha_crf,
                A_e_lj_l,  B_e_lj_l, A_e_crf_l, B_e_crf_l,
                A_de_lj_l, B_de_lj_l, A_de_crf_l, B_de_crf_l,
                lam);

            DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
            DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  storage.energies.A_lj_energy.size() " << conf.current().energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(it->i) << " energy group2 " 
                   << topo.atom_energy_group(it->j));

            conf.current().energies.A_lj_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_e_lj_l;

            conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_e_crf_l;

            conf.current().perturbed_energy_derivatives.A_lj_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_de_lj_l;

            conf.current().perturbed_energy_derivatives.A_crf_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_de_crf_l;
            DEBUG(8, "\ndone with storing energies ");
          } //all 101 lambda points done
        } // done with extended TI
      // ANITA
      
            break;
          }
          case simulation::pol_lj_crf_func : {
            /*double A_f6, A_f12;
            m_perturbed_nonbonded_term.
	      pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
					  A_lj->c6, A_lj->c12,
					  B_lj->c6, B_lj->c12,
					  A_qi, B_qi, A_qj, B_qj,
					  topo.coscharge(it->i),
					  topo.coscharge(it->j),
					  alpha_lj, alpha_crf,
					  A_f_pol, A_f6, A_f12,
					  A_e_lj, A_e_crf, A_de_lj, A_de_crf);
            // now combine everything
            A_f_pol_vec(0) = (A_f_pol[0] + A_f6 + A_f12) * r;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*double A_f6, A_f12;
            m_perturbed_nonbonded_term.
              pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                                          A_lj->c6, A_lj->c12,
                                          B_lj->c6, B_lj->c12,
                                          A_qi, B_qi, A_qj, B_qj,
                                          topo.coscharge(it->i),
                                          topo.coscharge(it->j),
                                          alpha_lj, alpha_crf,
                                          A_f_pol, A_f6, A_f12,
                                          A_e_lj, A_e_crf, A_de_lj, A_de_crf);
            // now combine everything
            A_f_pol_vec(0) = A_f_pol[0] * rm + (A_f6 + A_f12) * r;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      }
      else{ // not perturbed
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double A_f1 = 0.0;
            DEBUG(7, "non-perturbed interaction");
            m_nonbonded_term.
            lj_crf_interaction(r, A_lj->c6, A_lj->c12,
			       A_q, A_f1, A_e_lj,
			       A_e_crf);
            
            //A_de_lj = A_de_crf = 0.0;
            A_f = A_f1 * r * m_perturbed_nonbonded_term.A_lj_lambda_n();

	    A_de_lj = - topo.lambda_exp() * m_perturbed_nonbonded_term.A_lj_lambda_n_1() * A_e_lj; //again: calculate before scaling A_e_lj with lambda!
	    A_de_crf = - topo.lambda_exp() * m_perturbed_nonbonded_term.A_crf_lambda_n_1() * A_e_crf;
	    
	    
	      
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
          
          conf.current().energies.A_lj_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += A_e_lj; 
          conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += A_e_crf;

        }
      } // ANITA 

	    A_e_lj *= m_perturbed_nonbonded_term.A_lj_lambda_n();
	    A_e_crf *= m_perturbed_nonbonded_term.A_crf_lambda_n();

            break;
          }
          case simulation::pol_lj_crf_func : {
           /* m_nonbonded_term.
            pol_lj_crf_interaction(r, rp1, rp2, rpp, A_lj->c6, A_lj->c12,
                    topo.charge(it->i), topo.charge(it->j),
                    topo.coscharge(it->i), topo.coscharge(it->j),
                    A_f_pol, A_e_lj, A_e_crf);
            A_f_pol_vec(0) = A_f_pol[0] * r;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*m_nonbonded_term.
            pol_off_lj_crf_interaction(r, rm,  rp1, rp2, rpp, A_lj->c6, A_lj->c12,
                    topo.charge(it->i), topo.charge(it->j),
                    topo.coscharge(it->i), topo.coscharge(it->j),
                    A_f_pol, A_e_lj, A_e_crf);
            A_f_pol_vec(0) = A_f_pol[0] * r;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      }
 
      DEBUG(7, "\tcalculated interaction state A:\n\t\tf: "
	    << math::v2s(A_f) << " e_lj: " 
	    << A_e_lj << " e_crf: " << A_e_crf 
	    << " de_lj: " << A_de_lj 
	    << " de_crf: " << A_de_crf);
      
      break;
    case 2: // 1,4 interaction
      DEBUG(7, "\tlj-parameter state A cs6=" << A_lj->cs6 
	    << " cs12=" << A_lj->cs12);
      DEBUG(7, "\tcharges state A i*j = " << A_q);
      
      if (is_perturbed){
        DEBUG(7, "perturbed 1,4 interaction");
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double A_f1 = 0.0, A_f6 = 0.0, A_f12 = 0.0;
            
            m_perturbed_nonbonded_term.
	      lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
				      0, 0,
				      A_q, 0,
				      alpha_lj, alpha_crf,
				      A_f1, A_f6, A_f12, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
            A_f = (A_f1 + A_f6 + A_f12) * r;
            
      // ANITA
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
//        if ( sim.param().precalclam.nr_lambdas ) { 
          double A_e_lj_l = 0.0, B_e_lj_l = 0.0, A_e_crf_l = 0.0, B_e_crf_l = 0.0,
              A_de_lj_l = 0.0, B_de_lj_l = 0.0, A_de_crf_l = 0.0, B_de_crf_l = 0.0;

          // determine lambda stepsize from min,max and nr of lambdas
          double lambda_step = (sim.param().precalclam.max_lam - 
                   sim.param().precalclam.min_lam) / 
                   (sim.param().precalclam.nr_lambdas-1);

          //loop over nr_lambdas
          for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){ 

            // determine current lambda for this index
            double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

            // start the calculations
            m_perturbed_nonbonded_term.
            lj_crf_soft_interaction_ext(r, A_lj->cs6, A_lj->cs12,
                0, 0, A_q, 0, alpha_lj, alpha_crf,
                A_e_lj_l,  B_e_lj_l, A_e_crf_l, B_e_crf_l,
                A_de_lj_l, B_de_lj_l, A_de_crf_l, B_de_crf_l,
                lam);

            DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
            DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  storage.energies.A_lj_energy.size() " << conf.current().energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(it->i) << " energy group2 " 
                   << topo.atom_energy_group(it->j));

            conf.current().energies.A_lj_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_e_lj_l;

            conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_e_crf_l;

            conf.current().perturbed_energy_derivatives.A_lj_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_de_lj_l;

            conf.current().perturbed_energy_derivatives.A_crf_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += A_de_crf_l;
            DEBUG(8, "\ndone with storing energies ");
          } //all 101 lambda points done
        } // done with extended TI
      // ANITA
            break;
          }
          case simulation::pol_lj_crf_func : {
            /*double A_f6, A_f12;
            m_perturbed_nonbonded_term.
            pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
					A_lj->cs6, A_lj->cs12,
					B_lj->cs6, B_lj->cs12,
					A_qi, B_qi, A_qj, B_qj,
					topo.coscharge(it->i),
					topo.coscharge(it->j),
					alpha_lj, alpha_crf,
					A_f_pol, A_f6, A_f12,
					A_e_lj, A_e_crf, A_de_lj, A_de_crf);
            

            A_f_pol_vec(0) = (A_f_pol[0] + A_f6 + A_f12) * r;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	     */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*double A_f6, A_f12;
            m_perturbed_nonbonded_term.
            pol_off_lj_crf_soft_interaction(r,rm, rp1, rp2, rpp,
                                        A_lj->cs6, A_lj->cs12,
                                        B_lj->cs6, B_lj->cs12,
                                        A_qi, B_qi, A_qj, B_qj,
                                        topo.coscharge(it->i),
                                        topo.coscharge(it->j),
                                        alpha_lj, alpha_crf,
                                        A_f_pol, A_f6, A_f12,
                                        A_e_lj, A_e_crf, A_de_lj, A_de_crf);


            A_f_pol_vec(0) = A_f_pol[0] * rm + (A_f6 + A_f12) * r;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	     */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      } else { // non perturbed
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double A_f1 = 0.0;
            DEBUG(7, "non-perturbed 1,4 interaction");
            m_nonbonded_term.
            lj_crf_interaction(r, A_lj->cs6, A_lj->cs12,
			       A_q, A_f1, A_e_lj,
			       A_e_crf);
            //A_de_lj = A_de_crf = 0.0;

            A_f = A_f1 * r * m_perturbed_nonbonded_term.A_lj_lambda_n();

	    A_de_lj = - topo.lambda_exp() * m_perturbed_nonbonded_term.A_lj_lambda_n_1() * A_e_lj;
	    A_de_crf = - topo.lambda_exp() * m_perturbed_nonbonded_term.A_crf_lambda_n_1() * A_e_crf;
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
          
          conf.current().energies.A_lj_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += A_e_lj;
          conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += A_e_crf;
            
        }
      } // ANITA 

	    A_e_lj *= m_perturbed_nonbonded_term.A_lj_lambda_n();
	    A_e_crf *= m_perturbed_nonbonded_term.A_crf_lambda_n();


            break;
          }
          case simulation::pol_lj_crf_func : {
           /* m_nonbonded_term.
            pol_lj_crf_interaction(r, rp1, rp2, rpp, A_lj->cs6, A_lj->cs12,
				   topo.charge(it->i), topo.charge(it->j),
				   topo.coscharge(it->i), topo.coscharge(it->j),
				   A_f_pol, A_e_lj, A_e_crf);
            A_f_pol_vec(0) = A_f_pol[0] * r;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*m_nonbonded_term.
            pol_off_lj_crf_interaction(r, rm,  rp1, rp2, rpp, A_lj->cs6, A_lj->cs12,
                                   topo.charge(it->i), topo.charge(it->j),
                                   topo.coscharge(it->i), topo.coscharge(it->j),
                                   A_f_pol, A_e_lj, A_e_crf);
            A_f_pol_vec(0) = A_f_pol[0] * rm;
            A_f_pol_vec(1) = A_f_pol[1] * rp1;
            A_f_pol_vec(2) = A_f_pol[2] * rp2;
            A_f_pol_vec(3) = A_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      }
 
      DEBUG(7, "\tcalculated interaction state A:\n\t\tf: "
	    << math::v2s(A_f) << " e_lj: " 
	    << A_e_lj << " e_crf: " << A_e_crf 
	    << " de_lj: " << A_de_lj 
	    << " de_crf: " << A_de_crf);
      
      break;
      // --------------------------
  }
  
  // interaction in state B
  // ======================
  switch(it->B_type){
    // --------------------------
    case 0: // excluded
      B_e_lj = B_e_crf = B_de_lj = B_de_crf = 0.0;
      B_f = 0.0;
      DEBUG(7, "excluded in B");
      if(sim.param().nonbonded.rf_excluded){
        if (is_perturbed){
          switch(t_interaction_spec::interaction_func){
            case simulation::lj_crf_func : {
              m_perturbed_nonbonded_term.
              rf_soft_interaction(r, 0, B_q, alpha_crf,
				  B_f, B_e_crf, B_de_crf);
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
        double A_e_rf = 0.0, B_e_rf = 0.0, A_de_rf = 0.0, B_de_rf = 0.0;

        // determine lambda stepsize from min,max and nr of lambdas
        double lambda_step = (sim.param().precalclam.max_lam -
                 sim.param().precalclam.min_lam) /
                 (sim.param().precalclam.nr_lambdas-1);

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

          // determine current lambda for this index
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
          m_perturbed_nonbonded_term.rf_soft_interaction_ext(r, 0, B_q,
                  alpha_crf, A_e_rf, B_e_rf,
                  A_de_rf, B_de_rf, lam);
          conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += B_e_rf;
          conf.current().perturbed_energy_derivatives.B_crf_energy[lam_index]
            [topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += B_de_rf;

        }
      } // ANITA 

              break;
            }
            case simulation::pol_lj_crf_func : 
            case simulation::pol_off_lj_crf_func : {
              /*m_perturbed_nonbonded_term.
              pol_rf_soft_interaction(r, rp1, rp2, rpp,
				      A_qi, A_qj, B_qi, B_qj,
				      topo.coscharge(it->i), topo.coscharge(it->j),
				      alpha_crf, B_f_pol, B_e_crf, B_de_crf);
	     */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
              break;
            }
            default: io::messages.add("Perturbed_Nonbonded_Pair",
                    "interaction function not implemented",
                    io::message::critical);
          }
        } else { // non perturbed
          switch(t_interaction_spec::interaction_func){
            case simulation::lj_crf_func : {
              m_nonbonded_term.rf_interaction(r, B_q, B_f, B_e_crf);

	      B_f = B_f * m_perturbed_nonbonded_term.B_crf_lambda_n();

	      B_de_crf = topo.lambda_exp() * m_perturbed_nonbonded_term.B_crf_lambda_n_1() * B_e_crf;
	      
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

          conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += B_e_crf;

        }
      } // ANITA 


	      B_e_crf *=  m_perturbed_nonbonded_term.B_crf_lambda_n();


              break;
            }
            case simulation::pol_lj_crf_func : 
            case simulation::pol_off_lj_crf_func : {
              /*m_nonbonded_term.
              pol_rf_interaction(r, rp1, rp2, rpp, topo.charge()(it->i),
				 topo.charge()(it->j), topo.coscharge(it->i),
				 topo.coscharge(it->j), B_f_pol_vec, B_e_crf);
	      */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
              break;
            }
            default: io::messages.add("Perturbed_Nonbonded_Pair",
                    "interaction function not implemented",
                    io::message::critical);
          }
        }
      }

      break;
      // --------------------------
    case 1: // normal interaction
      DEBUG(7, "\tlj-parameter state B c6=" << B_lj->c6 
	    << " c12=" << B_lj->c12);
      DEBUG(7, "\tcharges state B i*j = " << B_q);
      
      if (is_perturbed){
	DEBUG(7, "perturbed interaction");
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double B_f1 = 0.0, B_f6 = 0.0, B_f12 = 0.0;
            
            m_perturbed_nonbonded_term.
            lj_crf_soft_interaction(r, 0, 0,
			            B_lj->c6, B_lj->c12,				    
				    0, B_q,
				    alpha_lj, alpha_crf,
				    B_f1, B_f6, B_f12,
				    B_e_lj, B_e_crf, B_de_lj, B_de_crf);
            
            B_f = (B_f1 + B_f6 + B_f12) * r;
           // ANITA
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
//        if ( sim.param().precalclam.nr_lambdas ) { 
          double A_e_lj_l = 0.0, B_e_lj_l = 0.0, A_e_crf_l = 0.0, B_e_crf_l = 0.0,
              A_de_lj_l = 0.0, B_de_lj_l = 0.0, A_de_crf_l = 0.0, B_de_crf_l = 0.0;

          // determine lambda stepsize from min,max and nr of lambdas
          double lambda_step = (sim.param().precalclam.max_lam -
                   sim.param().precalclam.min_lam) /
                   (sim.param().precalclam.nr_lambdas-1);

          //loop over nr_lambdas
          for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

            // determine current lambda for this index
            double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

            // start the calculations
            m_perturbed_nonbonded_term.
            lj_crf_soft_interaction_ext(r, 0, 0, B_lj->c6, B_lj->c12,
                0, B_q, alpha_lj, alpha_crf,
                A_e_lj_l,  B_e_lj_l, A_e_crf_l, B_e_crf_l,
                A_de_lj_l, B_de_lj_l, A_de_crf_l, B_de_crf_l,
                lam);

            DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
            DEBUG(8, "\n  B_e_lj " << B_e_lj << "\n  lambda index " << lam_index <<
                   "\n  conf.current().energies.B_lj_energy.size() " << conf.current().energies.B_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(it->i) << " energy group2 "
                   << topo.atom_energy_group(it->j));

            conf.current().energies.B_lj_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_e_lj_l;

            conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_e_crf_l;

            conf.current().perturbed_energy_derivatives.B_lj_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_de_lj_l;

            conf.current().perturbed_energy_derivatives.B_crf_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_de_crf_l;
            DEBUG(8, "\ndone with storing energies ");
          } //all 101 lambda points done
        } // done with extended TI
      // ANITA

            break;
          }
          case simulation::pol_lj_crf_func : {
            /*double B_f6, B_f12;
            m_perturbed_nonbonded_term.
            pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
					A_lj->c6, A_lj->c12,
					B_lj->c6, B_lj->c12,
					A_qi, B_qi, A_qj, B_qj,
					topo.coscharge(it->i),
					topo.coscharge(it->j),
					alpha_lj, alpha_crf,
					B_f_pol, B_f6, B_f12,
					B_e_lj, B_e_crf, B_de_lj, B_de_crf);
            // now combine everything
            B_f_pol_vec(0) = (B_f_pol[0] + B_f6 + B_f12) * r;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*double B_f6, B_f12;
            m_perturbed_nonbonded_term.
            pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                                        A_lj->c6, A_lj->c12,
                                        B_lj->c6, B_lj->c12,
                                        A_qi, B_qi, A_qj, B_qj,
                                        topo.coscharge(it->i),
                                        topo.coscharge(it->j),
                                        alpha_lj, alpha_crf,
                                        B_f_pol, B_f6, B_f12,
                                        B_e_lj, B_e_crf, B_de_lj, B_de_crf);
            // now combine everything
            B_f_pol_vec(0) = B_f_pol[0] * rm + (B_f6 + B_f12) * r;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
	     */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      }
      else{ // not perturbed
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double B_f1 = 0.0;
            DEBUG(7, "non-perturbed interaction");
            m_nonbonded_term.
            lj_crf_interaction(r, B_lj->c6, B_lj->c12,
			       B_q, B_f1, B_e_lj,
			       B_e_crf);
            
	  //  B_de_lj = B_de_crf = 0.0;

            B_f = B_f1 * r * m_perturbed_nonbonded_term.B_lj_lambda_n();

	    B_de_lj = topo.lambda_exp() * m_perturbed_nonbonded_term.B_lj_lambda_n_1() * B_e_lj;
	    B_de_crf = topo.lambda_exp() * m_perturbed_nonbonded_term.B_crf_lambda_n_1() * B_e_crf;
	    
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){


          conf.current().energies.B_lj_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += B_e_lj; 
          conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += B_e_crf; 

        }
      } // ANITA 

	    B_e_lj *= m_perturbed_nonbonded_term.B_lj_lambda_n();
	    B_e_crf *=  m_perturbed_nonbonded_term.B_crf_lambda_n();


            break;
          }
          case simulation::pol_lj_crf_func : {
            /*m_nonbonded_term.
            pol_lj_crf_interaction(r, rp1, rp2, rpp, B_lj->c6, B_lj->c12,
				   topo.charge(it->i), topo.charge(it->j),
				   topo.coscharge(it->i), topo.coscharge(it->j),
				   B_f_pol, B_e_lj, B_e_crf);
            B_f_pol_vec(0) = B_f_pol[0] * r;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*m_nonbonded_term.
            pol_off_lj_crf_interaction(r, rm,  rp1, rp2, rpp, B_lj->c6, B_lj->c12,
                                   topo.charge(it->i), topo.charge(it->j),
                                   topo.coscharge(it->i), topo.coscharge(it->j),
                                   B_f_pol, B_e_lj, B_e_crf);
            B_f_pol_vec(0) = B_f_pol[0] * r;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      }
 
      DEBUG(7, "\tcalculated interaction state B:\n\t\tf: "
	    << math::v2s(B_f) << " e_lj: " 
	    << B_e_lj << " e_crf: " << B_e_crf 
	    << " de_lj: " << B_de_lj 
	    << " de_crf: " << B_de_crf);
      
      break;
    case 2: // 1,4 interaction
      DEBUG(7, "\tlj-parameter state B cs6=" << B_lj->cs6 
	    << " cs12=" << B_lj->cs12);
      DEBUG(7, "\tcharges state B i*j = " << B_q);
      
      if (is_perturbed){
        DEBUG(7, "perturbed 1,4 interaction");
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double B_f1 = 0.0, B_f6 = 0.0, B_f12 = 0.0;
            
            m_perturbed_nonbonded_term.
            lj_crf_soft_interaction(r, 0, 0,
				    B_lj->cs6, B_lj->cs12,
				    0, B_q,
				    alpha_lj, alpha_crf,
				    B_f1, B_f6, B_f12, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
            B_f = (B_f1 + B_f6 + B_f12) * r;
            
          // ANITA
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){ 
          double A_e_lj_l = 0.0, B_e_lj_l = 0.0, A_e_crf_l = 0.0, B_e_crf_l = 0.0,
              A_de_lj_l = 0.0, B_de_lj_l = 0.0, A_de_crf_l = 0.0, B_de_crf_l = 0.0;

          // determine lambda stepsize from min,max and nr of lambdas
          double lambda_step = (sim.param().precalclam.max_lam -
                   sim.param().precalclam.min_lam) /
                   (sim.param().precalclam.nr_lambdas-1);

          //loop over nr_lambdas
          for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

            // determine current lambda for this index
            double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

            // start the calculations
            m_perturbed_nonbonded_term.
            lj_crf_soft_interaction_ext(r, 0, 0, B_lj->cs6, B_lj->cs12,
                 0, B_q, alpha_lj, alpha_crf,
                A_e_lj_l,  B_e_lj_l, A_e_crf_l, B_e_crf_l,
                A_de_lj_l, B_de_lj_l, A_de_crf_l, B_de_crf_l,
                lam);

            DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
            DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  conf.current().energies.A_lj_energy.size() " << conf.current().energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(it->i) << " energy group2 "
                   << topo.atom_energy_group(it->j));

            conf.current().energies.B_lj_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_e_lj_l;

            conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_e_crf_l;

            conf.current().perturbed_energy_derivatives.B_lj_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_de_lj_l;

            conf.current().perturbed_energy_derivatives.B_crf_energy
                    [lam_index][topo.atom_energy_group(it->i)]
                    [topo.atom_energy_group(it->j)] += B_de_crf_l;
            DEBUG(8, "\ndone with storing energies ");
          } //all 101 lambda points done
        } // done with extended TI
      // ANITA

            break;
          }
          case simulation::pol_lj_crf_func : {
            /*double B_f6, B_f12;
            m_perturbed_nonbonded_term.
            pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
					A_lj->cs6, A_lj->cs12,
					B_lj->cs6, B_lj->cs12,
					A_qi, B_qi, A_qj, B_qj,
					topo.coscharge(it->i),
					topo.coscharge(it->j),
					alpha_lj, alpha_crf,
					B_f_pol, B_f6, B_f12,
					B_e_lj, B_e_crf, B_de_lj, B_de_crf);
            
            //--------------------------------------------------
            // interactions have been calculated
            //--------------------------------------------------
            
            // now combine everything
            B_f_pol_vec(0) = (B_f_pol[0] + B_f6 + B_f12) * r;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*double B_f6, B_f12;
            m_perturbed_nonbonded_term.
            pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                                        A_lj->cs6, A_lj->cs12,
                                        B_lj->cs6, B_lj->cs12,
                                        A_qi, B_qi, A_qj, B_qj,
                                        topo.coscharge(it->i),
                                        topo.coscharge(it->j),
                                        alpha_lj, alpha_crf,
                                        B_f_pol, B_f6, B_f12,
                                        B_e_lj, B_e_crf, B_de_lj, B_de_crf);

            //--------------------------------------------------
            // interactions have been calculated
            //--------------------------------------------------

            // now combine everything
            B_f_pol_vec(0) = B_f_pol[0] * rm + (B_f6 + B_f12) * r;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
   	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      } else { // non perturbed
        switch(t_interaction_spec::interaction_func){
          case simulation::lj_crf_func : {
            double B_f1 = 0.0;
            DEBUG(7, "non-perturbed 1,4 interaction");
            m_nonbonded_term.
	      lj_crf_interaction(r, B_lj->cs6, B_lj->cs12,
				 B_q, B_f1, B_e_lj,
				 B_e_crf);
            B_de_lj = B_de_crf = 0.0;

            B_f = B_f1 * r * m_perturbed_nonbonded_term.B_lj_lambda_n();

	    B_de_lj = topo.lambda_exp() * m_perturbed_nonbonded_term.B_lj_lambda_n_1() * B_e_lj;
	    B_de_crf = topo.lambda_exp() * m_perturbed_nonbonded_term.B_crf_lambda_n_1() * B_e_crf;
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

          conf.current().energies.B_lj_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += B_e_lj;
          conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(it->i)]
            [topo.atom_energy_group(it->j)] += B_e_crf;
        }
      } // ANITA 

	    B_e_lj *= m_perturbed_nonbonded_term.B_lj_lambda_n();
	    B_e_crf *= m_perturbed_nonbonded_term.B_crf_lambda_n();


            break;
          }
          case simulation::pol_lj_crf_func : {
            /*m_nonbonded_term.
	      pol_lj_crf_interaction(r, rp1, rp2, rpp, B_lj->cs6, B_lj->cs12,
				     topo.charge(it->i), topo.charge(it->j),
				     topo.coscharge(it->i), topo.coscharge(it->j),
				     B_f_pol, B_e_lj, B_e_crf);
            B_f_pol_vec(0) = B_f_pol[0] * r;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
	    */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
            break;
          }
          case simulation::pol_off_lj_crf_func : {
            /*m_nonbonded_term.
              pol_off_lj_crf_interaction(r, rm, rp1, rp2, rpp, B_lj->cs6, B_lj->cs12,
                                     topo.charge(it->i), topo.charge(it->j),
                                     topo.coscharge(it->i), topo.coscharge(it->j),
                                     B_f_pol, B_e_lj, B_e_crf);
            B_f_pol_vec(0) = B_f_pol[0] * rm;
            B_f_pol_vec(1) = B_f_pol[1] * rp1;
            B_f_pol_vec(2) = B_f_pol[2] * rp2;
            B_f_pol_vec(3) = B_f_pol[3] * rpp;
	    */
            break;
          }
          default: io::messages.add("Perturbed_Nonbonded_Pair",
                  "interaction function not implemented",
                  io::message::critical);
        }
      }
 
      DEBUG(7, "\tcalculated interaction state B:\n\t\tf: "
	    << math::v2s(B_f) << " e_lj: " 
	    << B_e_lj << " e_crf: " << B_e_crf 
	    << " de_lj: " << B_de_lj 
	    << " de_crf: " << B_de_crf);
      
      break;
      // --------------------------
  }

  
  // A/B done
  
  if (perturbation_details::do_scaling){
    if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func) {
      io::messages.add("Perturbed_Nonbonded_Pair",
		       "interaction function not implemented",
		       io::message::critical);
      return;
    }
    if (t_interaction_spec::interaction_func == simulation::pol_off_lj_crf_func) {
      io::messages.add("Perturbed_Nonbonded_Pair",
                       "interaction function not implemented",
                       io::message::critical);
      return;
    }
    
    // check whether we need to do scaling
    // based on energy groups
    DEBUG(7, "scaled interaction: (perturbed pair) " << it->i << " - " << it->j);
    std::pair<int, int> energy_group_pair(topo.atom_energy_group(it->i),
					  topo.atom_energy_group(it->j));

    if (topo.energy_group_scaling().count(energy_group_pair)){
      
      // YES, we do scale the interactions!

      std::pair<double, double> scaling_pair =
	topo.energy_group_scaling()[energy_group_pair];
      A_f      *= scaling_pair.first;
      A_e_lj   *= scaling_pair.first;
      A_e_crf  *= scaling_pair.first;
      A_de_lj  *= scaling_pair.first;
      A_de_crf *= scaling_pair.first;
      B_f      *= scaling_pair.second;
      B_e_lj   *= scaling_pair.second;
      B_e_crf  *= scaling_pair.second;
      B_de_lj  *= scaling_pair.second;
      B_de_crf *= scaling_pair.second;
    }
  } // end of scaling

  // now combine everything
  DEBUG(10, "B_lj_l: " << m_perturbed_nonbonded_term.B_lj_lambda() <<
	" B_lj_ln: " << m_perturbed_nonbonded_term.B_lj_lambda_n() <<
	" A_lj_l: " << m_perturbed_nonbonded_term.A_lj_lambda() <<
	" A_lj_ln: " << m_perturbed_nonbonded_term.A_lj_lambda_n());
  /*
  if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func) {
     switch(it->A_type) {
      case 0 : // excluded
        if (is_perturbed) {
          // pol rf soft interaction: format already fine. do nothing
        } else {
          A_f = A_f_pol_vec(0) + A_f_pol_vec(1) + A_f_pol_vec(2) + A_f_pol_vec(3);
        }
        break;
      case 1 : // normal
      case 2 : // one four
        if (is_perturbed) {
          A_f = A_f_pol_vec(0) + A_f_pol_vec(1) + A_f_pol_vec(2) + A_f_pol_vec(3);
        } else {
          A_f = A_f_pol[0]*r + A_f_pol[1]*rp1 + A_f_pol[2]*rp2 + A_f_pol[3]*rpp;
        }
        break;
    }    
    
    switch(it->B_type) {
      case 0 : // excluded
        if (is_perturbed) {
          // pol rf soft interaction: format already fine. do nothing
        } else {
          B_f = B_f_pol_vec(0) + B_f_pol_vec(1) + B_f_pol_vec(2) + B_f_pol_vec(3);
        }
        break;
      case 1 : // normal
      case 2 : // one four
        if (is_perturbed) {
          B_f = B_f_pol_vec(0) + B_f_pol_vec(1) + B_f_pol_vec(2) + B_f_pol_vec(3);
        } else {
          B_f = B_f_pol[0]*r + B_f_pol[1]*rp1 + B_f_pol[2]*rp2 + B_f_pol[3]*rpp;
        }
        break;
    }   
  }
   */

 
/* OLD: lambdas were used twice for perturbed interactions, 
  	since e&f are already scaled by lambda in the perturbed_nonbonded functions. 
	non-perturbed-nonbonded interactions, however, are not scaled in their functions.


  // now combine state A and B
  // for the forces we take A_lj_lambda_n and B_lj_lambda_n, because we checked that this is the
  // same as A_crf_lambda_n and B_crf_lambda_n before

 if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func
      || t_interaction_spec::interaction_func == simulation::pol_off_lj_crf_func) {
    for(unsigned int i = 0; i < 4; ++i) {
      f_pol_vec(i) = m_perturbed_nonbonded_term.B_lj_lambda_n() * B_f_pol_vec(i) +
	m_perturbed_nonbonded_term.A_lj_lambda_n() * A_f_pol_vec(i);
    }
    f = f_pol_vec(0) +  f_pol_vec(1) +  f_pol_vec(2) +  f_pol_vec(3);
  } else {
    f = m_perturbed_nonbonded_term.B_lj_lambda_n() * B_f +
      m_perturbed_nonbonded_term.A_lj_lambda_n() * A_f;
  }
  
  e_lj   = m_perturbed_nonbonded_term.B_lj_lambda_n() * B_e_lj   + 
    m_perturbed_nonbonded_term.A_lj_lambda_n() * A_e_lj;
  e_crf  = m_perturbed_nonbonded_term.B_crf_lambda_n() * B_e_crf  + 
    m_perturbed_nonbonded_term.A_crf_lambda_n() * A_e_crf;

  de_lj  = m_perturbed_nonbonded_term.B_crf_lambda_n() * B_de_lj  + 
    m_perturbed_nonbonded_term.A_crf_lambda_n() * A_de_lj  
    + topo.lambda_exp() * m_perturbed_nonbonded_term.B_lj_lambda_n_1() * B_e_lj  
    - topo.lambda_exp() * m_perturbed_nonbonded_term.A_lj_lambda_n_1() * A_e_lj;
  
  de_crf = m_perturbed_nonbonded_term.B_crf_lambda_n() * B_de_crf + 
    m_perturbed_nonbonded_term.A_crf_lambda_n() * A_de_crf 
    + topo.lambda_exp() * m_perturbed_nonbonded_term.B_crf_lambda_n_1() * B_e_crf 
    - topo.lambda_exp() * m_perturbed_nonbonded_term.A_crf_lambda_n_1() * A_e_crf;
*/

//   NEW: nonbonded-non-perturbed interactions are scaled in-line,
//   and state A and B are added up here:

//forces:
  if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func
      || t_interaction_spec::interaction_func == simulation::pol_off_lj_crf_func) {
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);
  } else {
    f = B_f + A_f;
  }
//energies:
  e_lj   = B_e_lj + A_e_lj;
  e_crf  = B_e_crf + A_e_crf;
  
  de_lj  = B_de_lj  + A_de_lj;
  de_crf = B_de_crf + A_de_crf;


  conf.current().force(it->i) += f;
  conf.current().force(it->j) -= f;
  
  if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func
     || t_interaction_spec::interaction_func == simulation::pol_off_lj_crf_func) {
   /* for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
        conf.current().virial_tensor(a, b) += r(a) *( f_pol_vec(0)(b) +
                                               f_pol_vec(1)(b) + 
                                               f_pol_vec(2)(b) +
                                               f_pol_vec(3)(b));
  */
	      io::messages.add("Perturbed_Nonbonded_Pair", 
		      "polarization: interaction function not implemented",
                       io::message::critical);

  } else {
    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
        conf.current().virial_tensor(a, b) += r(a) * f(b);
  }

  
  DEBUG(7, "\tatomic virial done");

  DEBUG(7, "A_lj_lnm: " << m_perturbed_nonbonded_term.A_lj_lambda_n_1()
	<< " B_lj_lnm: " << m_perturbed_nonbonded_term.B_lj_lambda_n_1());
  DEBUG(7, "\tcalculated interaction:\n\t\tf: " << math::v2s(f) << " e_lj: " 
	<< e_lj << " e_crf: " << e_crf << " de_lj: " << de_lj 
	<< " de_crf: " << de_crf);
  
  // energy
  //assert(m_storage.energies.lj_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.energies.lj_energy.size() >
  //   topo.atom_energy_group(j));
  //assert(m_storage.energies.crf_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.energies.crf_energy.size() >
  //   topo.atom_energy_group(j));
  
  conf.current().energies.lj_energy
    [topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)] += e_lj;
  
  conf.current().energies.crf_energy
    [topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)] += e_crf;
  
  DEBUG(7, "\tenergy i and j " << topo.atom_energy_group(it->i)
	<< " " << topo.atom_energy_group(it->j));
  
  //assert(m_storage.perturbed_energy_derivatives.lj_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.perturbed_energy_derivatives.lj_energy.size() >
  //   topo.atom_energy_group(j));
  //assert(m_storage.perturbed_energy_derivatives.crf_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.perturbed_energy_derivatives.crf_energy.size() >
  //   topo.atom_energy_group(j));
  
  conf.current().perturbed_energy_derivatives.
    lj_energy[topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)]
    += de_lj;
  
  conf.current().perturbed_energy_derivatives.
    crf_energy[topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)]
    += de_crf;
  
}

