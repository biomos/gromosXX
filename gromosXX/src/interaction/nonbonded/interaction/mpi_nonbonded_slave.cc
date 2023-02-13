/**
 * @file mpi_nonbonded_slave.cc
 * methods of MPI_Nonbonded_Slave
 */
#ifdef XXMPI
#include <mpi.h>
#endif

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../simulation/parameter.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction.h"
#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"

#include "../../../interaction/nonbonded/interaction/storage.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_set.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_term.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_pair.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_set.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_set.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_interaction.h"
#include "../../../interaction/nonbonded/interaction/mpi_nonbonded_slave.h"

#include "../../../util/prepare_virial.h"

#include "../../../util/debug.h"
#include "../../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::MPI_Nonbonded_Slave::MPI_Nonbonded_Slave(Pairlist_Algorithm *pa)
  : Nonbonded_Interaction(pa)
{
}

/**
 * Destructor.
 * @bug change destruction of nonbonded set to be standard - conform!
 */
interaction::MPI_Nonbonded_Slave::~MPI_Nonbonded_Slave()
{
  DEBUG(7, "MPI_Nonbonded_Slave::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::MPI_Nonbonded_Slave::calculate_interactions
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(4, "MPI_Nonbonded_Slave::calculate_interactions");

  m_timer.start(sim);

#ifdef XXMPI

  // check if we want to calculate nonbonded
  // might not be necessary if multiple time-stepping is enabled
  
  int steps = sim.param().multistep.steps;
  if (steps == 0) steps = 1;
  
  if ((sim.steps() % steps) == 0){
    // std::cout << "MULTISTEP: full calculation\n";

    int rank = sim.mpiControl().threadID;
    int num_threads = sim.mpiControl().numberOfThreads;
    
    // distribute the positions
    assert(((double *) &conf.current().pos(conf.current().pos.size()-1)(0)) -
	   ((double *) &conf.current().pos(0)(0)) == int((conf.current().pos.size() - 1)*3));
    
    // std::cerr << "slave: receiving pos" << std::endl;
    MPI_Bcast(&conf.current().pos(0)(0),
			  conf.current().pos.size() * 3, 
			  MPI::DOUBLE,
			   sim.mpiControl().masterID, sim.mpiControl().comm);
   
    // receiving charges for QM/MM
    MPI_Bcast(&topo.charge()[0],
    topo.num_atoms(),
        MPI::DOUBLE,
        sim.mpiControl().masterID, sim.mpiControl().comm);

	


    // std::cerr << "slave: receiving box" << std::endl;
    MPI_Bcast(&conf.current().box(0)(0),
	                  9,
                          MPI::DOUBLE,
                           sim.mpiControl().masterID, sim.mpiControl().comm);
    // std::cerr << "slave: clfldsa; d" << std::endl;
    
    // bcast lambda for slow growth and chemical monte carlo
    MPI_Bcast(&topo.lambda(),
                          1,
                          MPI::DOUBLE,
                           sim.mpiControl().masterID, sim.mpiControl().comm);
    // and also update the individual lambdas on the slave!
    topo.update_for_lambda();
 
    // do this on the master and on the slaves...
    m_pairlist_algorithm->prepare(topo, conf, sim);
    
    // prepare for the virial
    util::prepare_virial(topo, conf, sim);
    
    // calculate interactions for our rank
    DEBUG(8, " SLAVE::calculating nonbonded interactions (thread " 
	  << rank << " of " << num_threads << ")");
    
    m_nonbonded_set[0]->calculate_interactions(topo, conf, sim);
    
    // collect the forces, energies, energy-derivatives, virial
    // MPI::IN_PLACE ???
    MPI_Reduce(&m_nonbonded_set[0]->storage().force(0)(0),
			   NULL,
			   m_nonbonded_set[0]->storage().force.size() * 3,
			   MPI::DOUBLE,
			   MPI::SUM,
			    sim.mpiControl().masterID, sim.mpiControl().comm);
    const unsigned int ljs = conf.current().energies.lj_energy.size();
    if (sim.param().force.force_groups) {
      for (unsigned int i = 0; i < ljs; ++i) {
        for (unsigned int j = 0; j < ljs; ++j) {
          MPI_Reduce(&m_nonbonded_set[0]->storage().force_groups[i][j](0)(0),
                  NULL, m_nonbonded_set[0]->storage().force_groups[i][j].size() * 3,
                  MPI::DOUBLE, MPI::SUM,
                   sim.mpiControl().masterID, sim.mpiControl().comm);
        }
      }
    }
    std::vector<double> lj_scratch(ljs*ljs);
    std::vector<double> crf_scratch(ljs*ljs);
    std::vector<double> ls_real_scratch(ljs*ljs);
    
    for(unsigned int i = 0; i < ljs; ++i){
      for(unsigned int j = 0; j < ljs; ++j){
	lj_scratch[i*ljs + j] = 
	  m_nonbonded_set[0]->storage().energies.lj_energy[i][j];
	crf_scratch[i*ljs + j] = 
	  m_nonbonded_set[0]->storage().energies.crf_energy[i][j];
        ls_real_scratch[i*ljs + j] =
          m_nonbonded_set[0]->storage().energies.ls_real_energy[i][j];      
      }
    }
    MPI_Reduce(&lj_scratch[0],
			   NULL,
			   ljs * ljs,
			   MPI::DOUBLE,
			   MPI::SUM,
			    sim.mpiControl().masterID, sim.mpiControl().comm);
    MPI_Reduce(&crf_scratch[0],
			   NULL,
			   ljs * ljs,
			   MPI::DOUBLE,
			   MPI::SUM,
			    sim.mpiControl().masterID, sim.mpiControl().comm);
    MPI_Reduce(&ls_real_scratch[0],
			   NULL,
			   ljs * ljs,
			   MPI::DOUBLE,
			   MPI::SUM,
			    sim.mpiControl().masterID, sim.mpiControl().comm);
    if (sim.param().pcouple.virial){
      double * dvt2 = &m_nonbonded_set[0]->storage().virial_tensor(0,0);
      MPI_Reduce(dvt2,
			     NULL,
			     9,
			     MPI::DOUBLE,
			     MPI::SUM,
			      sim.mpiControl().masterID, sim.mpiControl().comm);
    }
    
    if (sim.param().perturbation.perturbation){
      
      for(unsigned int i = 0; i < ljs; ++i){
	for(unsigned int j = 0; j < ljs; ++j){
	  lj_scratch[i*ljs + j] = 
	    m_nonbonded_set[0]->storage().perturbed_energy_derivatives.lj_energy[i][j];
	  crf_scratch[i*ljs + j] = 
	    m_nonbonded_set[0]->storage().perturbed_energy_derivatives.crf_energy[i][j];
	}
      }
      MPI_Reduce(&lj_scratch[0],
			     NULL,
			     ljs * ljs,
			     MPI::DOUBLE,
			     MPI::SUM,
			      sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&crf_scratch[0],
			     NULL,
			     ljs * ljs,
			     MPI::DOUBLE,
			     MPI::SUM,
			      sim.mpiControl().masterID, sim.mpiControl().comm);
      
    }

    // ANITA
    unsigned int nr_lambdas = sim.param().precalclam.nr_lambdas;
    for(unsigned int l = 0; l < nr_lambdas; ++l){
      std::vector<double> A_lj_scratch(ljs*ljs);
      std::vector<double> B_lj_scratch(ljs*ljs);
      std::vector<double> A_crf_scratch(ljs*ljs);
      std::vector<double> B_crf_scratch(ljs*ljs);

      std::vector<double> A_dlj_scratch(ljs*ljs);
      std::vector<double> B_dlj_scratch(ljs*ljs);
      std::vector<double> A_dcrf_scratch(ljs*ljs);
      std::vector<double> B_dcrf_scratch(ljs*ljs);

      for(unsigned int i = 0; i < ljs; ++i){
        for(unsigned int j = 0; j < ljs; ++j){
          A_lj_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().energies.A_lj_energy[l][i][j];
          B_lj_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().energies.B_lj_energy[l][i][j];
          A_crf_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().energies.A_crf_energy[l][i][j];
          B_crf_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().energies.B_crf_energy[l][i][j];

          A_dlj_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().perturbed_energy_derivatives.A_lj_energy[l][i][j];
          B_dlj_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().perturbed_energy_derivatives.B_lj_energy[l][i][j];
          A_dcrf_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().perturbed_energy_derivatives.A_crf_energy[l][i][j];
          B_dcrf_scratch[i*ljs + j] =
             m_nonbonded_set[0]->storage().perturbed_energy_derivatives.B_crf_energy[l][i][j];
        }
      }

      MPI_Reduce(&A_lj_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&B_lj_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&A_crf_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&B_crf_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);

      MPI_Reduce(&A_dlj_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&B_dlj_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&A_dcrf_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&B_dcrf_scratch[0],
              NULL,
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);

    } // ANITA
    
    if (sim.param().pairlist.print &&
	 (!(sim.steps() % sim.param().pairlist.skip_step))){
      
      DEBUG(7, "print pairlist...");
      std::ofstream os("slave.pl", std::ios::app);
      os << "rank " << rank << " of " << num_threads << std::endl;
      print_pairlist(topo, conf, sim, os);
    }
    
    if (sim.param().eds.eds){
      const unsigned int numstates = sim.param().eds.numstates;
      
      assert(m_nonbonded_set[0]->storage().virial_tensor_endstates.size() == numstates);
      assert(m_nonbonded_set[0]->storage().force_endstates.size() == numstates);
      assert(m_nonbonded_set[0]->storage().energies.eds_vi.size() == numstates);
      assert(m_nonbonded_set[0]->storage().energies.eds_vi_shift_extra_orig.size() == numstates);
      assert(m_nonbonded_set[0]->storage().energies.eds_vi_shift_extra_phys.size() == numstates);
      
      // reduce energies of endstates
      MPI_Reduce(&m_nonbonded_set[0]->storage().energies.eds_vi[0],
              NULL,
              numstates,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);

      // reduce extra energies of endstates
      MPI_Reduce(&m_nonbonded_set[0]->storage().energies.eds_vi_shift_extra_orig[0],
              NULL,
              numstates,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);

      // reduce extra energies of endstates
      MPI_Reduce(&m_nonbonded_set[0]->storage().energies.eds_vi_shift_extra_phys[0],
              NULL,
              numstates,
              MPI::DOUBLE,
              MPI::SUM,
               sim.mpiControl().masterID, sim.mpiControl().comm);
      
      // reduce virial tensors of endstates
      if (sim.param().pcouple.virial){
        
        double * dvt2 = &m_nonbonded_set[0]->storage().virial_tensor_endstates[0](0, 0);
        MPI_Reduce(dvt2,
                NULL,
                9 * numstates,
                MPI::DOUBLE,
                MPI::SUM,
                 sim.mpiControl().masterID, sim.mpiControl().comm);
        
      }
      
      for(unsigned int state = 0; state < numstates; state++){
        
        // reduce forces of endstates
        MPI_Reduce(&m_nonbonded_set[0]->storage().force_endstates[state](0)(0),
                NULL,
                m_nonbonded_set[0]->storage().force_endstates[state].size() * 3,
                MPI::DOUBLE,
                MPI::SUM,
                 sim.mpiControl().masterID, sim.mpiControl().comm);
         /*
        // reduce energies of endstates
        MPI::COMM_WORLD.Reduce(&m_nonbonded_set[0]->storage().energies.eds_vi[state],
                NULL,
                1,
                MPI::DOUBLE,
                MPI::SUM,
                0);
          
        // reduce virial tensors of endstates
        if (sim.param().pcouple.virial){
          
          double * dvt2 = &m_nonbonded_set[0]->storage().virial_tensor_endstates[state](0, 0);
          MPI::COMM_WORLD.Reduce(dvt2,
                  NULL,
                  9,
                  MPI::DOUBLE,
                  MPI::SUM,
                  0);
           
        }
         */
      } // loop over states
    } // eds
    
    if (sim.param().nonbonded.method == simulation::el_p3m ||
        sim.param().nonbonded.method == simulation::el_ewald) {
      MPI_Reduce(&m_nonbonded_set[0]->storage().energies.ls_kspace_total,
            NULL, 1, MPI::DOUBLE, MPI::SUM,  sim.mpiControl().masterID, sim.mpiControl().comm);
    }

    if(sim.param().force.interaction_function == simulation::lj_shifted_crf_corr_func){
      std::vector<double> shift_extra_orig_scratch(ljs*ljs);
      std::vector<double> shift_extra_phys_scratch(ljs*ljs);
      
      for(unsigned int i = 0; i < ljs; ++i){
        for(unsigned int j = 0; j < ljs; ++j){
          shift_extra_orig_scratch[i*ljs + j] = 
            m_nonbonded_set[0]->storage().energies.shift_extra_orig[i][j];
          shift_extra_phys_scratch[i*ljs + j] = 
            m_nonbonded_set[0]->storage().energies.shift_extra_phys[i][j];
        }
      }
      MPI_Reduce(&shift_extra_orig_scratch[0],
           NULL,
           ljs * ljs,
           MPI::DOUBLE,
           MPI::SUM,
           sim.mpiControl().masterID, sim.mpiControl().comm);

      MPI_Reduce(&shift_extra_phys_scratch[0],
           NULL,
           ljs * ljs,
           MPI::DOUBLE,
           MPI::SUM,
           sim.mpiControl().masterID, sim.mpiControl().comm);
    }
      
    ////////////////////////////////////////////////////
    // end of multiple time stepping: calculate
    ////////////////////////////////////////////////////
  }
  else{
    // std::cout << "MULTISTEP: no recalculation...\n";
  }

  DEBUG(6, "MPI_Nonbonded_Slave::calculate_interactions done");
  m_timer.stop();
  
  return 0;
#else
  std::cerr << "using MPI code without MPI defined..." << std::endl;
  return E_ILLEGAL;
#endif
}

/**
 * initialize the arrays
 * need to override to pass MPI rank, size to nonbonded set
 */
int interaction::MPI_Nonbonded_Slave::init
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::ostream & os,
 bool quiet
 )
{
#ifdef XXMPI
  int rank =  sim.mpiControl().threadID;
  int num_threads = sim.mpiControl().numberOfThreads;
  
  // initialise the pairlist...
  m_pairlist_algorithm->init(topo, conf, sim, os, quiet);
  
  if (sim.param().nonbonded.method != simulation::el_reaction_field) {
    conf.lattice_sum().init(topo, sim);
  }  

  DEBUG(15, "MPI_Nonbonded_Slave::initialize");
  m_nonbonded_set.clear();

  if (sim.param().perturbation.perturbation){
    
    // only one set per MPI process
    m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
							  m_parameter, rank, num_threads));
  }
  else if (sim.param().eds.eds){
    m_nonbonded_set.push_back(new Eds_Nonbonded_Set(*m_pairlist_algorithm,
            m_parameter, rank, num_threads));
  }
  else{
    // only one set per MPI process
    m_nonbonded_set.push_back(new Nonbonded_Set(*m_pairlist_algorithm, 
						m_parameter, rank, num_threads));
  }
  
  if (sim.param().multicell.multicell)
    m_nonbonded_set[0]->init(topo.multicell_topo(), conf, sim, os, quiet);
  else
    m_nonbonded_set[0]->init(topo, conf, sim, os, quiet);
  
  if (check_special_loop(topo, conf, sim, os, quiet) != 0) {
    io::messages.add("special solvent loop check failed", "Nonbonded_Interaction",
            io::message::error);
  }
  return 0;
  
#else
  std::cerr << "MPI: MPI_Nonbonded_Master::init but MPI not defined" << std::endl;
  return 1;
#endif

}


//***************************************************************************
// helper functions 
//***************************************************************************


