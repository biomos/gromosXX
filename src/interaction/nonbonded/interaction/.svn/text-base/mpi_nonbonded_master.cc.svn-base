/**
 * @file mpi_nonbonded_master.cc
 * methods of MPI_Nonbonded_Master
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
#include "../../../interaction/nonbonded/interaction/mpi_nonbonded_master.h"

#include "../../../util/debug.h"
#include "../../../util/error.h"



#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::MPI_Nonbonded_Master::MPI_Nonbonded_Master(Pairlist_Algorithm *pa)
: Nonbonded_Interaction(pa)
{
}

/**
 * Destructor.
 * @bug change destruction of nonbonded set to be standard - conform!
 */
interaction::MPI_Nonbonded_Master::~MPI_Nonbonded_Master() 
{
  DEBUG(7, "MPI_Nonbonded_Master::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::MPI_Nonbonded_Master::
calculate_interactions(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) 
{
  DEBUG(4, "MPI_Nonbonded_Master::calculate_interactions");
  
  if (sim.param().multicell.multicell){
    io::messages.add("MPI code with multiple unit cell simulations not implemented",
            "mpi_nonbonded_interaction",
            io::message::critical);
    return 1;
  }
  
  m_timer.start();
  
#ifdef XXMPI
  
  // check if we want to calculate nonbonded
  // might not be necessary if multiple time-stepping is enabled
  
  int steps = sim.param().multistep.steps;
  if (steps == 0) steps = 1;
  
  if ((sim.steps() % steps) == 0){
    // std::cout << "MULTISTEP: full calculation\n";
    
    int rank = MPI::COMM_WORLD.Get_rank();
    int num_threads = MPI::COMM_WORLD.Get_size();
    
    // --------------------------------------------------
    // distribute the positions
    
    // if this is true, the positions (and forces) can be copied
    // directly from theVArray
    assert(((double *) &conf.current().pos(conf.current().pos.size()-1)(0)) -
            ((double *) &conf.current().pos(0)(0)) == int((conf.current().pos.size() - 1)*3));
    
    // std::cerr << "master: bcast pos" << std::endl;
    MPI::COMM_WORLD.Bcast(&conf.current().pos(0)(0),
            conf.current().pos.size() * 3,
            MPI::DOUBLE,
            0);
   
	 // bcast charges (for QM/MM)
	                    MPI::COMM_WORLD.Bcast(&topo.charge()[0],
	                                           topo.num_atoms(),
	                                              MPI::DOUBLE,
	                                         0);



    // don't forget the box (or are you stupid or what????)
    // std::cerr << "master: bcast box" << std::endl;
    MPI::COMM_WORLD.Bcast(&conf.current().box(0)(0),
            9,
            MPI::DOUBLE,
            0);
    
    // bcast lambda for slow growth and chemical monte carlo
    MPI::COMM_WORLD.Bcast(&topo.lambda(),
            1,
            MPI::DOUBLE,
            0);
    
    // std::cerr << "ready to calc" << std::endl;
    
    // --------------------------------------------------
    // calculate interactions
    
    DEBUG(8, "calculating nonbonded interactions (thread "
            << rank << " of " << num_threads << ")");
    
    // do this on the master and on the slaves...
    m_pairlist_algorithm->prepare(topo, conf, sim);
    
    m_nonbonded_set[0]->calculate_interactions(topo, conf, sim);
    
    // collect the forces, energies, energy-derivatives, virial
    // MPI::IN_PLACE ???
    m_storage.force = m_nonbonded_set[0]->storage().force;
    MPI::COMM_WORLD.Reduce(&m_storage.force(0)(0),
            &m_nonbonded_set[0]->storage().force(0),
            m_nonbonded_set[0]->storage().force.size() * 3,
            MPI::DOUBLE,
            MPI::SUM,
            0);
    

    const unsigned int ljs = conf.current().energies.lj_energy.size();
    if (sim.param().force.force_groups) {
      for (unsigned int i = 0; i < ljs; ++i) {
        for (unsigned int j = 0; j < ljs; ++j) {
          m_storage.force_groups[i][j] = m_nonbonded_set[0]->storage().force_groups[i][j];
          MPI::COMM_WORLD.Reduce(&m_storage.force_groups[i][j](0)(0),
            &m_nonbonded_set[0]->storage().force_groups[i][j](0),
            m_nonbonded_set[0]->storage().force_groups[i][j].size() * 3,
            MPI::DOUBLE,
            MPI::SUM,
            0);
        }
      }
    }
    std::vector<double> lj_scratch(ljs*ljs), rlj_scratch(ljs*ljs);
    std::vector<double> crf_scratch(ljs*ljs), rcrf_scratch(ljs*ljs);
    std::vector<double> ls_real_scratch(ljs*ljs), rls_real_scratch(ljs*ljs);
    
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
    MPI::COMM_WORLD.Reduce(&lj_scratch[0],
            &rlj_scratch[0],
            ljs * ljs,
            MPI::DOUBLE,
            MPI::SUM,
            0);
    MPI::COMM_WORLD.Reduce(&crf_scratch[0],
            &rcrf_scratch[0],
            ljs * ljs,
            MPI::DOUBLE,
            MPI::SUM,
            0);
    MPI::COMM_WORLD.Reduce(&ls_real_scratch[0],
            &rls_real_scratch[0],
            ljs * ljs,
            MPI::DOUBLE,
            MPI::SUM,
            0);
    for(unsigned int i = 0; i < ljs; ++i){
      for(unsigned int j = 0; j < ljs; ++j){
        m_nonbonded_set[0]->storage().energies.lj_energy[i][j] = rlj_scratch[i*ljs + j];
        m_nonbonded_set[0]->storage().energies.crf_energy[i][j] = rcrf_scratch[i*ljs + j];
        m_nonbonded_set[0]->storage().energies.ls_real_energy[i][j] = rls_real_scratch[i*ljs + j];
      }
    }
    
    if (sim.param().pcouple.virial){
      
      m_storage.virial_tensor = m_nonbonded_set[0]->storage().virial_tensor;
      double * dvt = &m_storage.virial_tensor(0,0);
      double * dvt2 = &m_nonbonded_set[0]->storage().virial_tensor(0,0);
      MPI::COMM_WORLD.Reduce(dvt,
              dvt2,
              9,
              MPI::DOUBLE,
              MPI::SUM,
              0);
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
      MPI::COMM_WORLD.Reduce(&lj_scratch[0],
              &rlj_scratch[0],
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
              0);
      MPI::COMM_WORLD.Reduce(&crf_scratch[0],
              &rcrf_scratch[0],
              ljs * ljs,
              MPI::DOUBLE,
              MPI::SUM,
              0);
      
      for(unsigned int i = 0; i < ljs; ++i){
        for(unsigned int j = 0; j < ljs; ++j){
          m_nonbonded_set[0]->storage().perturbed_energy_derivatives.lj_energy[i][j]
                  = rlj_scratch[i*ljs + j];
          m_nonbonded_set[0]->storage().perturbed_energy_derivatives.crf_energy[i][j]
                  = rcrf_scratch[i*ljs + j];
        }
      }
    }
    
    if (sim.param().eds.eds){
      const unsigned int numstates = sim.param().eds.numstates;
      
      assert(m_storage.virial_tensor_endstates.size() == numstates);
      assert(m_nonbonded_set[0]->storage().virial_tensor_endstates.size() == numstates);
      assert(m_storage.force_endstates.size() == numstates);
      assert(m_nonbonded_set[0]->storage().force_endstates.size() == numstates);
      assert(m_storage.energies.eds_vi.size() == numstates);
      assert(m_nonbonded_set[0]->storage().energies.eds_vi.size() == numstates);
      
      m_storage.energies.eds_vi  = m_nonbonded_set[0]->storage().energies.eds_vi;
      m_storage.force_endstates = m_nonbonded_set[0]->storage().force_endstates;
      m_storage.virial_tensor_endstates = m_nonbonded_set[0]->storage().virial_tensor_endstates;
       
      // reduce energies of endstates
      MPI::COMM_WORLD.Reduce(&m_storage.energies.eds_vi[0],
              &m_nonbonded_set[0]->storage().energies.eds_vi[0],
              numstates,
              MPI::DOUBLE,
              MPI::SUM,
              0);
      
      // reduce virial tensors of endstates
      if (sim.param().pcouple.virial){
        double * dvt = &m_storage.virial_tensor_endstates[0](0, 0);
        double * dvt2 = &m_nonbonded_set[0]->storage().virial_tensor_endstates[0](0, 0);
        MPI::COMM_WORLD.Reduce(dvt,
                dvt2,
                9 * numstates,
                MPI::DOUBLE,
                MPI::SUM,
                0);
      }

      for(unsigned int state = 0; state < numstates; state++){
        
        // reduce forces of endstates
        MPI::COMM_WORLD.Reduce(&m_storage.force_endstates[state](0)(0),
                &m_nonbonded_set[0]->storage().force_endstates[state](0),
                m_nonbonded_set[0]->storage().force_endstates[state].size() * 3,
                MPI::DOUBLE,
                MPI::SUM,
                0);
         /*
        // reduce energies of endstates
        MPI::COMM_WORLD.Reduce(&m_storage.energies.eds_vi[state],
                &m_nonbonded_set[0]->storage().energies.eds_vi[state],
                1,
                MPI::DOUBLE,
                MPI::SUM,
                0);
                 
        // reduce virial tensors of endstates
        if (sim.param().pcouple.virial){
          double * dvt = &m_storage.virial_tensor_endstates[state](0, 0);
          double * dvt2 = &m_nonbonded_set[0]->storage().virial_tensor_endstates[state](0, 0);
          MPI::COMM_WORLD.Reduce(dvt,
                  dvt2,
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
      double sum = 0.0;
      MPI::COMM_WORLD.Reduce(&m_nonbonded_set[0]->storage().energies.ls_kspace_total,
            &sum, 1, MPI::DOUBLE, MPI::SUM, 0);
      m_nonbonded_set[0]->storage().energies.ls_kspace_total = sum;
    }
     
    ////////////////////////////////////////////////////
    // end of multiple time stepping: calculate
    ////////////////////////////////////////////////////
    
    if (sim.param().pairlist.print &&
    (!(sim.steps() % sim.param().pairlist.skip_step))){
      
      DEBUG(7, "print pairlist...");
      std::ofstream os("server.pl", std::ios::app);
      os << "rank " << rank << " of " << num_threads << std::endl;
      print_pairlist(topo, conf, sim, os);
    }
    
  }
  else{
    // std::cout << "MULTISTEP: no recalculation...\n";
  }
  
  DEBUG(6, "sets are done, adding things up...");
  store_set_data(topo, conf, sim);
  
  
  DEBUG(6, "Nonbonded_Interaction::calculate_interactions done");
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
int interaction::MPI_Nonbonded_Master::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) 
{
#ifdef XXMPI
  int rank = MPI::COMM_WORLD.Get_rank();
  int num_threads = MPI::COMM_WORLD.Get_size();
  
  // initialise the pairlist...
  m_pairlist_algorithm->init(topo, conf, sim, os, quiet);
  
  if (sim.param().nonbonded.method != simulation::el_reaction_field) {
    conf.lattice_sum().init(topo, sim);
  }

  
  DEBUG(15, "MPI_Nonbonded_Master::initialize");
  m_nonbonded_set.clear();
  
  if (sim.param().perturbation.perturbation){
    
    // only one set per MPI process
    m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
            m_parameter, rank, num_threads));
  }
  else if (sim.param().eds.eds){
    m_nonbonded_set.push_back(new Eds_Nonbonded_Set(*m_pairlist_algorithm,
            m_parameter, rank, num_threads));
    
    m_storage.force_endstates.resize(sim.param().eds.numstates);
    m_storage.virial_tensor_endstates.resize(sim.param().eds.numstates);
    m_storage.energies.eds_vi.resize(sim.param().eds.numstates);
    
    for(unsigned int i = 0; i < m_storage.force_endstates.size(); i++){
      m_storage.force_endstates[i].resize(topo.num_atoms());
    }
     
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
  
  m_storage.force.resize(topo.num_atoms());
  m_storage.energies.
  resize(unsigned(conf.current().energies.bond_energy.size()),
          unsigned(conf.current().energies.kinetic_energy.size()));
  
  if (sim.param().force.force_groups) {
    m_storage.force_groups.resize(unsigned(conf.current().energies.bond_energy.size()),
            std::vector<math::VArray>(unsigned(conf.current().energies.bond_energy.size()), 
            math::VArray(topo.num_atoms(), math::Vec(0.0, 0.0, 0.0))));
  }
  
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


