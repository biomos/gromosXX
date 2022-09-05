/**
 * @file qmmm_interaction.cc
 * Implements QMMM interaction
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../interaction/interaction.h"

#include "../../../util/debug.h"
#include "../../../util/error.h"

#include "../../../math/boundary_checks.h"

#include "mm_atom.h"
#include "qm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "nonbonded/qmmm_nonbonded_set.h"
#include "qmmm_interaction.h"

/*
#ifdef OMP
#include <omp.h>
#endif
*/
#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::QMMM_Interaction * interaction::QMMM_Interaction::qmmm_ptr = nullptr;

interaction::QMMM_Interaction::QMMM_Interaction() : Interaction("QMMM Interaction")
                                                  , m_parameter()
                                                  , m_set_size(1)
                                                  , m_rank(0)
                                                  , m_size(1)
                                                  , m_worker(nullptr)
                                                  , m_qm_zone(nullptr)
                                                  , m_qm_buffer(nullptr)
  {
#ifdef XXMPI
    m_rank = MPI::COMM_WORLD.Get_rank();
    m_size = MPI::COMM_WORLD.Get_size();
#endif
  qmmm_ptr = this;
}

interaction::QMMM_Interaction::~QMMM_Interaction() {
  for (unsigned int i = 0; i < m_qmmm_nonbonded_set.size(); ++i) {
    DEBUG(12, "deleting set " << i);
    delete m_qmmm_nonbonded_set[i];
    m_qmmm_nonbonded_set[i] = nullptr;
  }
  if (m_worker != nullptr)
    DEBUG(12, "deleting QM Worker");
    delete m_worker;
  if (m_qm_zone != nullptr)
    DEBUG(12, "deleting QM Zone");
    delete m_qm_zone;
  if (m_qm_buffer != nullptr)
    DEBUG(12, "deleting QM Buffer");
    delete m_qm_buffer;
}

/**
 * This should be separate class AddRemove : public QM_Link
 * 
configuration::Configuration interaction::QMMM_Interaction::AddRemove(topology::Topology &topo,
                                                                      configuration::Configuration &conf,
                                                                      simulation::Simulation &sim) {
    // create a conf for capped system. Required for the evaluation of the forces on the external
    //point charges in the MOPAC worker
    //
    configuration::Configuration qmmm_conf = conf;
    qmmm_conf.current().force=0.0;
    math::Vec posCap,dR;
    unsigned int m1,q1;
    const double rch=0.109; //Carbon-Hydrogen bond length
    const bool verbose = false;
    unsigned int pi=0;
        for (std::vector< std::pair<unsigned int,unsigned int> >::const_iterator
                     it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it,++pi  )
        {
            q1=it->first;
            m1=it->second;
            dR = conf.current().pos(m1) - conf.current().pos(q1);
            posCap = (rch/abs(dR)    )    * dR + conf.current().pos(q1) ;
            qmmm_conf.current().pos(m1) = posCap;
        }
    return qmmm_conf;
}
*/

/**
 * This should be separate class AddRemove2 : public QM_Link
 *
int interaction::QMMM_Interaction::AddRemove2(topology::Topology &topo,
                                                                      configuration::Configuration &conf,
                                                                      simulation::Simulation &sim) {
    math::Vec posCap,dQM_MM,unit_MM_QM;
    math::Vec FQM,FMM;
    double d_quot=0;
    double d_cart=0;
    unsigned int m1,q1;
    const double rch=0.109; //Carbon-Hydrogen bond length
    const bool verbose = false;
    unsigned int pi=0;
    for (std::vector< std::pair<unsigned int,unsigned int> >::const_iterator
                 it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it,++pi  )
    {
        q1=it->first;
        m1=it->second;

       // std::cout << "size qmpair" << topo.qm_mm_pair().size()<<  std::endl;
        dQM_MM = conf.current().pos(m1) - conf.current().pos(q1);
        unit_MM_QM=dQM_MM/abs(dQM_MM);
        posCap = (rch/abs(dQM_MM)    )    * dQM_MM + conf.current().pos(q1) ;
        d_quot=rch/abs(dQM_MM);
        //std::cout << "pi" << pi << "  abs(unit_mm)" << abs(unit_MM_QM) << std::endl;
        for (unsigned int i = 0; i < 3; i++) {
            d_cart=dQM_MM[i];
            FQM[i]=storage.force(q1)[i]+stor_link.force(pi)[i]*((1.-d_quot)+(d_cart*unit_MM_QM[i]*d_quot/abs(dQM_MM)));
            FMM[i]=stor_link.force(pi)[i]*((d_quot)-(d_cart*unit_MM_QM[i]*d_quot/abs(dQM_MM)));
            //FMM=conf.current().force(m1)[i]+stor_link.force(pi)[i]*((d_quot)-(d_cart*unit_MM_QM[i]*d_quot/abs(dQM_MM)));
        }
       // std::cout << "FQM " << stor_link.force(pi)[0] << "  " <<  stor_link.force(pi)[1] << "  "
        //          << stor_link.force(pi)[2] << std::endl;
        //  std::cout << "FQM " << FQM[0] << "  " << FQM[1] << "  " << FQM[2] << std::endl;
     //   std::cout << "FMM " << FMM[0] << "  " << FMM[1] << "  " << FMM[2] << std::endl;
 //       conf.current().force(q1)+=FQM;
//        conf.current().force(m1)+=FMM;
        storage.force(q1)=FQM; //otherwise doublecounting when QM-forces are added to conf.current().force
        storage.force(m1)=FMM; //otherwise doublecounting when QM-forces are added to conf.current().force
    }

}
// useless comment

*/
int interaction::QMMM_Interaction::scf_step(topology::Topology& topo,
                                            configuration::Configuration& conf,
                                            simulation::Simulation& sim) {
  // Update only COS in QM zone
  m_timer.start();
  DEBUG(15,"Updating COS in QM zone");
  m_timer.start("QM zone update");
  m_qm_zone->update_cos(topo, conf, sim);
  m_timer.stop("QM zone update");

  int err = 0;
  m_timer.start(m_worker->name());
  err = m_worker->run_QM(topo, conf, sim, *m_qm_zone);
  m_timer.stop(m_worker->name());
  if (err) return err;
  m_timer.stop();
  return 0;
}

void interaction::QMMM_Interaction::write_qm_data(topology::Topology& topo,
                                                  configuration::Configuration& conf,
                                                  const simulation::Simulation& sim) {
  DEBUG(15,"Writing QM data");
  m_timer.start("writing QM results");
  m_qm_zone->write(topo, conf, sim);
  m_timer.stop("writing QM results");
}


int interaction::QMMM_Interaction::calculate_interactions(topology::Topology& topo,
                                                          configuration::Configuration& conf,
                                                          simulation::Simulation& sim) {
  DEBUG(4, "QMMM_Interaction::calculate_interactions");
  int err = 0;

  // Do QMMM only on master (MPI)
  if (m_rank == 0) {
    m_timer.start();
    // Update QM Zone
    DEBUG(15,"Updating QM zone");
    m_timer.start("QM zone update");
    err = m_qm_zone->update(topo, conf, sim);
    m_timer.stop("QM zone update");
    if (err) return err;
    
    m_timer.start(m_worker->name());
    err = m_worker->run_QM(topo, conf, sim, *m_qm_zone);
    m_timer.stop(m_worker->name());
    if (err) return err;

    if (sim.param().qmmm.use_qm_buffer) {
      DEBUG(4, "Using QM buffer");
      /** If we are using buffer region, this region is treated as both QM and MM
       * We construct our system from 3 subsystems:
       * 1. Outer region (OR)
       * 2. Buffer region (BR)
       * 3. QM region (QR)
       * Now we treat the interactions as follows:
       * OR energies - classical MM interactions (possible addition of delta(BRQR - BR) from calculation in electrostatic embedding)
       * OR forces - classical MM interactions (possible addition of delta(BRQR - BR) from calculation in electrostatic embedding)
       * BR energies - classical MM interactions + delta(BRQR - BR)
       * BR forces - classical MM interactions + delta(BRQR - BR)
       * QR energies - delta(BRQR - BR) only
       * QR forces - delta(BRQR - BR) only
       * Energy term delta(BRQR - BR) from the twin QM calculation contains inseparable energy contributions
       * 
       * Now we can get delta(BRQR - BR) from
       * a) one evaluation of BRQR with NN trained on deltas
       * b) from 2 QM or NN evaluations and calculating the difference here
       * 
       */
      if (sim.param().qmmm.software != simulation::qm_nn
          || sim.param().qmmm.nn.model_type == simulation::nn_model_type_standard) {
        DEBUG(4, "Creating QM buffer for separate QM calculation");
        //create buffer zone for separate QM calculation and run it
        delete m_qm_buffer;
        m_qm_buffer = m_qm_zone->create_buffer_zone(topo, sim);
        m_timer.start(m_worker->name());
        err = m_worker->run_QM(topo, conf, sim, *m_qm_buffer);
        m_timer.stop(m_worker->name());
        if (err) return err;
        // Calculate QM energy and QM forces as difference
        /**
         * this->evaluate_buffered_qm(m_qm_zone, m_qm_buffer)
         */
        {

          DEBUG(7, "QM+Buffer Zone energy: " << m_qm_zone->QM_energy());
          DEBUG(7, "Buffer Zone energy: " << m_qm_buffer->QM_energy());
          DEBUG(7, "Delta: " << m_qm_zone->QM_energy() - m_qm_buffer->QM_energy());
          m_qm_zone->QM_energy() -= m_qm_buffer->QM_energy();

          for (std::set<interaction::QM_Atom>::const_iterator buf_it = m_qm_buffer->qm.begin();
                buf_it != m_qm_buffer->qm.end(); ++buf_it) {
            std::set<interaction::QM_Atom>::iterator qm_it = m_qm_zone->qm.find(*buf_it);
            DEBUG(10, "QM/Buffer atom " << qm_it->index << "/" << buf_it->index);
            DEBUG(10, "QM+Buffer force: " << math::v2s(qm_it->force));
            DEBUG(10, "Buffer force: " <<  math::v2s(buf_it->force));
            qm_it->force -= buf_it->force;
            DEBUG(10, "Delta: " <<  math::v2s(qm_it->force));
          }
          for (std::set<interaction::MM_Atom>::const_iterator buf_it = m_qm_buffer->mm.begin();
                buf_it != m_qm_buffer->mm.end(); ++buf_it) {
            std::set<interaction::MM_Atom>::iterator mm_it = m_qm_zone->mm.find(*buf_it);
            DEBUG(10, "QM/Buffer MM atom " << mm_it->index << "/" << buf_it->index);
            DEBUG(10, "QM+Buffer force: " <<  math::v2s(mm_it->force));
            DEBUG(10, "Buffer force: " <<  math::v2s(buf_it->force));
            mm_it->force -= buf_it->force;
            DEBUG(10, "Delta: " <<  math::v2s(mm_it->force));
          }
        }
      } else { // BuRNN model (a)
        DEBUG(1, "Skipping buffer zone calculation");
        DEBUG(1, "Expecting deltas directly from NN");
      }
    }
    
    if (sim.param().qmmm.qmmm != simulation::qmmm_polarisable) {
      // in polarisable embedding, we will write the data after electric field evaluation
      this->write_qm_data(topo, conf, sim);
    }

    if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical
        || sim.param().qmmm.qm_lj) {
      this->calculate_nonbonded(topo, conf, sim);
    }
    m_timer.stop();
  }
  return 0;
}

int interaction::QMMM_Interaction::calculate_nonbonded(topology::Topology& topo,
                                                       configuration::Configuration& conf,
                                                       simulation::Simulation& sim)
  {
  DEBUG(4, "QMMM_Interaction::calculate_nonbonded");

  // allow multistep - calculate nonbonded only every nth step
  //int steps = sim.param().multistep.steps;
  //if (steps == 0) steps = 1;

  //if ((sim.steps() % steps) == 0) {
    for (unsigned i = 0; i < m_set_size; ++i) {
      if(m_qmmm_nonbonded_set[i]->calculate_interactions(topo, conf, sim)) {
        m_timer.stop();
	      return 1;
      }
    }
  //}
  DEBUG(6, "sets are done, adding things up...");
  this->store_set_data(topo, conf, sim);
  
  if (sim.param().pairlist.print &&
      (!(sim.steps() % sim.param().pairlist.skip_step))) {
    DEBUG(7, "print QM-MM pairlist...");
    this->print_pairlist(topo);
  }

  DEBUG(6, "QMMM_Interaction::calculate_nonbonded done");
  return 0;
}

void interaction::QMMM_Interaction::get_electric_field(const simulation::Simulation& sim
                                                     , math::VArray & electric_field) {
  // Write electric field
  m_timer.start();
  m_timer.start("Electric field writing");
  m_qm_zone->electric_field(sim, electric_field);
  m_timer.stop("Electric field writing");
  m_timer.stop();
}

int interaction::QMMM_Interaction::init(topology::Topology& topo,
            configuration::Configuration& conf,
            simulation::Simulation& sim,
            std::ostream& os,
            bool quiet) {
  if (!quiet)
    os << "QMMM INTERACTION\n";
    // Initial checks
  if (sim.param().qmmm.cutoff > 0.0 && 
          !math::boundary_check_cutoff(conf.current().box, sim.param().boundary.boundary,
          sim.param().qmmm.cutoff)) {
    io::messages.add("The RCUTQ cutoff is too large for the size of the "
            "computational box.", "QMMM_Interaction", io::message::error);
    return 1;
  }

  if (m_rank == 0) {  
    // Create QM_Zone
    const int charge = sim.param().qmmm.qm_zone.charge + sim.param().qmmm.buffer_zone.charge;
    int sm = sim.param().qmmm.qm_zone.spin_mult;
    if (sim.param().qmmm.use_qm_buffer) {
      // Calculate combined spin multiplicity
      // number of unpaired spins of the QM zone
      const int spin_qm = sim.param().qmmm.qm_zone.spin_mult - 1;
      // number of unpaired spins of the buffer zone
      const int spin_buf = sim.param().qmmm.buffer_zone.spin_mult - 1;
      // consider no spin pairing between the QM and buffer zone
      sm = spin_qm + spin_buf + 1;
    }
    delete m_qm_zone;
    m_qm_zone = new interaction::QM_Zone(charge, sm);

    DEBUG(15,"QM Zone created");
    DEBUG(15,"Net charge: " << charge);
    DEBUG(15,"Spin multiplicity: " << sm);

    if (m_qm_zone->init(topo, conf, sim)) return 1;
    DEBUG(15,"QM Zone initialized");

    DEBUG(15,"Creating QM Worker");
    m_worker = interaction::QM_Worker::get_instance(sim);
    if (m_worker == nullptr || m_worker->init(topo, conf, sim, *(m_qm_zone))) {
      io::messages.add("Error initializing QM worker", "QMMM_Interaction", io::message::error);
      return 1;
    }
    DEBUG(15,"QM Worker initialized");
  }
  if (!quiet) {
    switch (sim.param().qmmm.qmmm) {
      case simulation::qmmm_mechanical:
        os << "\tmechanical";
        break;
      case simulation::qmmm_electrostatic:
        os << "\telectrostatic";
        break;
      case simulation::qmmm_polarisable:
        os << "\tpolarisable";
        break;
      default:
        os << "\tunknown";
    }
    os << " embedding scheme" << std::endl;
    if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical) {
      os << "\tcharges of QM atoms ";
      switch (sim.param().qmmm.qm_ch) {
        case simulation::qm_ch_constant:
          os << "constant from the topology";
          break;
        case simulation::qm_ch_dynamic:
          os << "updated every step from the QM calculation";
          break;
        default:
          os << "unknown";
          break;
      }
      os << std::endl;
    }
    os << "\tusing external ";
    switch (sim.param().qmmm.software) {
      case simulation::qm_mndo:
        os << "MNDO";
        break;
      case simulation::qm_dftb:
        os << "DFTB";
        break;
      case simulation::qm_mopac:
        os << "MOPAC";
        break;
      case simulation::qm_turbomole:
        os << "Turbomole";
        break;
      case simulation::qm_gaussian:
        os << "Gaussian";
        break;
      case simulation::qm_nn:
        os << "Schnet";
        break;
      case simulation::qm_orca:
        os << "Orca";
        break;
#ifdef WITH_XTB
      case simulation::qm_xtb:
        os << "XTB";
        break;
#endif
      default:
        os << "unknown";
        break;
    }
    os << " program package" << std::endl;
    os.precision(3);
    os << "\tunits conversion factors: ";
    m_worker->print_unit_factors(os);
    os << std::endl;
    if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical) {
      os << "\tQM-MM interactions will be treated classically" << std::endl
         << "\tusing standard cutoffs (RCUTP, RCUTL)" << std::endl
         << "\tRCUTQM value (" << sim.param().qmmm.cutoff << ") ignored";
    }
    else if (sim.param().qmmm.cutoff == 0.0) {
      os << "\tincluding all MM atoms in QM-MM interaction calculation";
    }
    else {
      if (sim.param().qmmm.atomic_cutoff) os << "\tatom-based cutoff";
      else os << "\tchargegroup-based cutoff";
    os << ": " << std::fabs(sim.param().qmmm.cutoff) << " nm";
    }
    os << std::endl;

    if (!topo.qmmm_link().empty()) {
      os << "\tusing link-atom scheme with capping atom" << std::endl
         << "\tdistance between QM link atom and capping atom: "
         << sim.param().qmmm.cap_length << std::endl
         << "\tQM-MM links:" << std::endl;
      for (std::set< std::pair<unsigned,unsigned> >::const_iterator
            it = topo.qmmm_link().begin()
          , to = topo.qmmm_link().end()
          ; it != to; ++it) {
        os << "\t\t" << it->first << "-" << it->second << std::endl;
      }
    }

    if (sim.param().qmmm.qm_lj)
      os << "\tLJ interactions between QM atoms enabled" << std::endl;
    else
      os << "\tLJ interactions between QM atoms disabled" << std::endl;
      
    if (!sim.param().qmmm.qm_constraint)
      os << "\tremoving QM-QM constraints" << std::endl;

    if (sim.param().qmmm.mm_scale > 0.0) {
      os << "\tMM point charges will be scaled using (2/pi)*atan(s*|R|) with s = " <<
        sim.param().qmmm.mm_scale << std::endl;
      os << "\t\t|R| is distance to the closest QM atom" << std::endl;
    }
    unsigned num_qm = 0;
    unsigned num_buffer = 0;
    for (unsigned i = 0; i < topo.num_atoms(); ++i) {
      if (topo.is_qm(i))
        ++num_qm;
      else if (topo.is_qm_buffer(i))
        ++num_buffer;
    }
    os << "\tQM zone: " << std::endl
       << "\t\tnet charge         : " << sim.param().qmmm.qm_zone.charge << std::endl
       << "\t\tspin multiplicity  : " << sim.param().qmmm.qm_zone.spin_mult << std::endl
       << "\t\tnumber of QM atoms : " << num_qm << std::endl;
    if (sim.param().qmmm.use_qm_buffer) {
      os << "\t";
      if (sim.param().qmmm.buffer_zone.cutoff)
        os << "adaptive ";
      else
        os << "static ";
      os << "buffer zone:" << std::endl;
      if (sim.param().qmmm.buffer_zone.cutoff)
        os << "\t\tcutoff                    : " << sim.param().qmmm.buffer_zone.cutoff <<std::endl;
      os <<   "\t\tnet charge                : " << sim.param().qmmm.buffer_zone.charge << std::endl
         <<   "\t\tspin multiplicity         : " << sim.param().qmmm.buffer_zone.spin_mult << std::endl
         <<   "\t\tnumber of buffer atoms    : ";
      if (sim.param().qmmm.buffer_zone.cutoff)
        os << "up to ";
      os << num_buffer << std::endl;
    }
    os <<     "\tnumber of QM-MM links\t      : " << m_qm_zone->link.size() << std::endl;

  }

  // Remove relevant bonded terms from topo
  DEBUG(15, "Removing bonded terms");
  this->remove_bonded_terms(topo, os, quiet);

  if (!quiet)
    os << "\n";

  // Remove exclusions containing QM-QM and QM-MM interactions
  DEBUG(15, "Removing QM-QM and QM-MM exclusions");
  this->modify_exclusions(topo, sim, os, quiet);

  // Create nonbonded set for LJ interactions
  if (sim.param().force.nonbonded_vdw
      && (sim.param().qmmm.qmmm > simulation::qmmm_mechanical // Will do LJ between QM and MM
        || sim.param().qmmm.qm_lj)) // Will do LJ between QM atoms
    {
    this->init_nonbonded(topo, conf, sim, os, quiet);
  }

  if (!quiet)
    os << "END" << std::endl;
  DEBUG(9, "QMMM init done");
  return 0;
}

int interaction::QMMM_Interaction::init_nonbonded(topology::Topology& topo,
                                                  configuration::Configuration& conf,
                                                  simulation::Simulation& sim,
                                                  std::ostream& os,
                                                  bool quiet) {
  DEBUG(9, "QMMM nonbonded init");
  if (!quiet)
    os << "\tQMMM nonbonded interaction:\n";

  if (sim.param().multicell.multicell) {
    io::messages.add("MULTICELL is not allowed with QMMM"
                    , "QMMM_Interaction", io::message::error);
    return 1;
  }

  m_qmmm_nonbonded_set.clear();

// OpenMP and MPI parallelization if necessary
/*
#ifdef OMP
  m_set_size *= omp_get_num_threads();
#endif
#ifdef XXMPI
  m_set_size *= MPI::COMM_WORLD.Get_size();*/
  for (unsigned i = 0; i < m_set_size; ++i) {
    m_qmmm_nonbonded_set.push_back(
          new QMMM_Nonbonded_Set(*(this->m_qm_zone), this->m_timer
                                , m_parameter, i, m_set_size));
  }

  if (!quiet)
    os << "\t\tcreated " << m_qmmm_nonbonded_set.size() << " set"
       <<  (m_qmmm_nonbonded_set.size() > 1 ? "s" : "") << std::endl;

  std::vector<QMMM_Nonbonded_Set *>::iterator
    it = m_qmmm_nonbonded_set.begin(),
    to = m_qmmm_nonbonded_set.end();

  bool q = quiet;
  for (; it != to; ++it) {
    (*it)->init(topo, conf, sim, os, q);
    // only print first time...
    q = true;
  }
  DEBUG(9, "QMMM nonbonded init done");
  return 0;
}

void interaction::QMMM_Interaction::remove_bonded_terms(
                                topology::Topology& topo
                              , std::ostream& os
                              , bool quiet)
  {
  // Remove bonded terms that will be defined by QM
  // Definitions to simplify the code
  typedef std::vector<topology::two_body_term_struct> twoBodyVec;
  typedef std::vector<topology::three_body_term_struct> threeBodyVec;
  typedef std::vector<topology::four_body_term_struct> fourBodyVec;
  typedef std::vector<topology::eight_body_term_struct> eightBodyVec;
  twoBodyVec& bonds = topo.solute().bonds();
  twoBodyVec& cgbonds = topo.solute().cgbonds();
  threeBodyVec& angles = topo.solute().angles();
  fourBodyVec& improper_dihedrals = topo.solute().improper_dihedrals();
  fourBodyVec& dihedrals = topo.solute().dihedrals();
  eightBodyVec& crossdihedrals = topo.solute().crossdihedrals();
  if (!quiet)
    os << "\tterms removed from topology:";
  // Counter for neat printing
  unsigned count = 0;
  if (!quiet && bonds.size())
    os << "\n\t\tbonds:\n";
  
  std::string link_err_msg = "across QM-MM boundary, but no QMMM link defined between atoms ";
  auto is_qm_or_buffer = [&topo](const unsigned i)-> bool {
    return (topo.is_qm(i) || topo.is_qm_buffer(i));
  };
  // Delete QM-QM bonded terms and check, if QM-MM links were properly defined
  for (twoBodyVec::iterator b_it = bonds.begin(); b_it != bonds.end(); ) {
    // If QM-QM or QM-Buf
    if ((topo.is_qm(b_it->i) && is_qm_or_buffer(b_it->j))
        || (is_qm_or_buffer(b_it->i) && topo.is_qm(b_it->j))
      ) {
      if (!quiet) {
        if (count == 0) os << "\t\t";
        os << (b_it->i + 1) << "-" << (b_it->j + 1) << " ";
        if (++count == 8) {
          os << "\n";
          count = 0;
        }
      }
      DEBUG(4,"Erased bond: " << b_it->i + 1 << " - " << b_it->j + 1);
      b_it = bonds.erase(b_it);
      continue;
    }
    // If Buf-Buf
    else if (topo.is_qm_buffer(b_it->i) && topo.is_qm_buffer(b_it->i))
      ;
    // If QM/Buf-MM or MM-QM/Buf and not in QM-MM link
    else if ((is_qm_or_buffer(b_it->i)
            && !topo.are_linked( b_it->i, b_it->j ))
          ||
            (is_qm_or_buffer(b_it->j)
            && !topo.are_linked( b_it->j, b_it->i ))
          )
      {
      std::ostringstream msg;
      msg << "Bonded interaction " << link_err_msg
          << b_it->i << " and " << b_it->j;
      io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
    }
    ++b_it;
  }

  if (!quiet && cgbonds.size())
    os << "\n\t\tcoarse-grained bonds:\n";
  // Delete coarse-grained bonds between QM-QM
  for (twoBodyVec::iterator b_it = cgbonds.begin(); b_it != cgbonds.end(); ) {
    // If QM-QM or QM-Buf
    if ((topo.is_qm(b_it->i) && is_qm_or_buffer(b_it->j))
        || (is_qm_or_buffer(b_it->i) && topo.is_qm(b_it->j))
      ) {
      if (!quiet)
        os << "\t\t" << (b_it->i + 1) << " " << (b_it->j + 1) << "\n";
      DEBUG(4,"Erased coarse-grained bond: " << b_it->i + 1 << "-" << b_it->j + 1);
      b_it = cgbonds.erase(b_it);
      continue;
    } else ++b_it;
  }

  if (!quiet && angles.size()) {
    os << "\n\t\tbond angles:\n";
    count = 0;
  }
  /* Delete (QM/Buf/MM)-QM-(QM/Buf/MM) bond-angle terms and check,
   * if QM-MM links were properly defined
   */
  for (threeBodyVec::iterator a_it = angles.begin(); a_it != angles.end(); ) {
    const unsigned i = a_it->i
                 , j = a_it->j
                 , k = a_it->k;
    const unsigned qm_count = topo.is_qm(i) + topo.is_qm(j) + topo.is_qm(k);
    if (qm_count == 0) { ++a_it; continue; }// MM-MM-MM
    if (topo.is_qm(j)) {
      a_it = angles.erase(a_it);
      if (!quiet) {
        if (count == 0) os << "\t\t";
        os << (i+1) << "-" << (j+1) << "-" << (k+1) << " ";
        if (++count == 5) {
          os << "\n";
          count = 0;
        }
      }
      DEBUG(4,"Erased bond-angle bending term: " << i+1 << "-" << j+1 << "-" << k+1);
      if (qm_count == 3) continue; // QM-QM-QM, nothing to check
    } else ++a_it;
    // Loop over neighbouring pairs and check, if QM-MM links were properly defined
    const std::vector<unsigned> indices = {i,j,k};
    for (std::vector<unsigned>::const_iterator
        i_it = indices.begin(), j_it = i_it + 1, i_to = indices.end() - 1
        ; i_it != i_to; ++i_it, ++j_it) {
      if ((   is_qm_or_buffer(*i_it)
          && !is_qm_or_buffer(*j_it)
          && !topo.are_linked(*i_it, *j_it))
        ||
          (  !is_qm_or_buffer(*i_it)
          &&  is_qm_or_buffer(*j_it)
          && !topo.are_linked(*j_it, *i_it)
        ))
        {
        std::ostringstream msg;
        msg << "Bond-angle bending interaction " << link_err_msg
          << *i_it + 1 << " and " << *j_it + 1;
        io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
      }
    }
  }

  if (!quiet && dihedrals.size()) {
    os << "\n\t\tdihedral angles:\n";
    count = 0;
  }
  // Delete (MM/Buf/QM)-QM-QM-(MM/Buf/QM) terms - dihedrals
  for (fourBodyVec::iterator d_it = dihedrals.begin(); d_it != dihedrals.end(); ) {
    const unsigned i = d_it->i
                 , j = d_it->j
                 , k = d_it->k
                 , l = d_it->l;
    const unsigned qm_count = topo.is_qm(i) + topo.is_qm(j) + topo.is_qm(k) + topo.is_qm(l);
    if (qm_count == 0) { ++d_it; continue; }
    if ((topo.is_qm(j) && topo.is_qm(k))) {
      d_it = dihedrals.erase(d_it);
      if (!quiet) {
        if (count == 0) os << "\t\t";
        os << (i+1) << "-" << (j+1) << "-" << (k+1) << "-" << (l+1) << " ";
        if (++count == 4) {
          os << "\n";
          count = 0;
        }
      }
      DEBUG(4,"Erased dihedral angle torsion term: "
            << i+1 << "-" << j+1 << "-" << k+1 << "-" << l+1);
      if (qm_count == 4) continue;
    } else ++d_it;
    // loop over neighbouring pairs and check, if QM-MM links were properly defined
    const std::vector<unsigned> indices = {i,j,k,l};
    for (std::vector<unsigned>::const_iterator
        i_it = indices.begin(), j_it = i_it + 1, to = indices.end()
        ; j_it != to; ++i_it, ++j_it) {
      if ((   is_qm_or_buffer(*i_it)
          && !is_qm_or_buffer(*j_it)
          && !topo.are_linked(*i_it, *j_it))
        ||
          (  !is_qm_or_buffer(*i_it)
          &&  is_qm_or_buffer(*j_it)
          && !topo.are_linked(*j_it, *i_it)
        ))
        {
        std::ostringstream msg;
        msg << "Dihedral angle torsion interaction " << link_err_msg
            << *i_it + 1 << " and " << *j_it + 1;
        io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
      }
    }
  }
  DEBUG(15,"Dihedral angle torsion terms done");
  
  if (!quiet && improper_dihedrals.size()) {
    os << "\n\t\timproper dihedral angles:\n";
    count = 0;
  }
  // Delete QM-(MM/Buf/QM)-(MM/Buf/QM)-(MM/Buf/QM) improper dihedral
  for (fourBodyVec::iterator
      d_it = improper_dihedrals.begin(), d_to; d_it != improper_dihedrals.end(); ) {
    const unsigned i = d_it->i
                 , j = d_it->j
                 , k = d_it->k
                 , l = d_it->l;
    const unsigned qm_count = topo.is_qm(i) + topo.is_qm(j) + topo.is_qm(k) + topo.is_qm(l);
    if (qm_count == 0) { ++d_it; continue; }
    if ( topo.is_qm(d_it->i) )
      {
      DEBUG(4,"Erased improper dihedral angle bending term: "
            << i + 1 << "-" << j + 1 << "-"
            << k + 1 << "-" << l + 1);
      d_it = improper_dihedrals.erase(d_it);
      if (!quiet) {
        if (count == 0) os << "\t\t";
        os << (i+1) << "-" << (j+1) << "-" << (k+1) << "-" << (l+1) << " ";
        if (++count == 4) {
          os << "\n";
          count = 0;
        }
      }
      if (qm_count == 4) continue;
    } else ++d_it;

    // loop over neighbouring pairs and check, if QM-MM links were properly defined
    const std::vector<unsigned> indices = {j,k,l};
    for (std::vector<unsigned>::const_iterator
        j_it = indices.begin(), to = indices.end()
        ; j_it != to; ++j_it) {
      if ((   is_qm_or_buffer(i)
          && !is_qm_or_buffer(*j_it)
          && !topo.are_linked(i, *j_it))
        ||
          (  !is_qm_or_buffer(i)
          &&  is_qm_or_buffer(*j_it)
          && !topo.are_linked(*j_it, i)
        ))
        {
        std::ostringstream msg;
        msg << "Improper dihedral angle bending interaction " << link_err_msg
            << i + 1 << " and " << *j_it + 1;
        io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
      }
    }
  }

  // Delete All-QM crossdihedrals
  for (eightBodyVec::iterator d_it = crossdihedrals.begin(); d_it != crossdihedrals.end(); ) {
    if (   (topo.is_qm(d_it->a))
        && (topo.is_qm(d_it->b))
        && (topo.is_qm(d_it->c))
        && (topo.is_qm(d_it->d))
        && (topo.is_qm(d_it->e))
        && (topo.is_qm(d_it->f))
        && (topo.is_qm(d_it->g))
        && (topo.is_qm(d_it->h))
      )
      {
      if (!quiet) {
      os << "\t\tremoved crossdihedral bending term: "
         << d_it->a + 1 << " " << d_it->b + 1 << " "
         << d_it->c + 1 << " " << d_it->d + 1 << " "
         << d_it->e + 1 << " " << d_it->f + 1 << " "
         << d_it->g + 1 << " " << d_it->h + 1 << std::endl;
      }
      d_it = crossdihedrals.erase(d_it);
    } else ++d_it;
  }
}

void interaction::QMMM_Interaction::modify_exclusions(
                                topology::Topology& topo
                              , const simulation::Simulation& sim
                              , std::ostream& os
                              , bool quiet)
    /*** MECHANICAL EMBEDDING */
    /** We need to remove:
     *    1. QM-QM one-four pairs - should be 0
     *    2. QM-QM exclusions - to avoid RF correction term
     *                        - they are excluded anyway within pairlist
     *                        - keep copy and calculate LJ, if QMLJ is on
     *    3. We need to do QM-QM LJ exceptions here (if QMLJ is ON)
     * * QM-MM pairs are done classically with full charges, exclusions and LJs
     */

    /** ELECTROSTATIC EMBEDDING */
    /** The same as mechanical embedding, but also remove and make a copy of:
     *    1. QM-MM one-four pairs - should be run here with LJ only
     *                            - this maybe makes sense only for
     *                              link-atom schemes
     *    2. QM-MM exclusions - to avoid RF and LS correction term
     *    3. QM-MM exceptions - should be done here
     */
  {
  DEBUG(4, "Removing exclusions of QM atoms");
  topo.qm_all_exclusion().clear();
  topo.qm_all_exclusion().resize(topo.num_atoms());
  // Remove or make a copy of QM-QM and QM-MM exclusions
  const bool use_qm_buffer = sim.param().qmmm.use_qm_buffer;
  
  /** (K)eep or (R)emove?
   *     QM  Buf  MM
   * QM   R    R  R/K
   * Buf  R    K   K
   * MM   R    K   K
   *
   * R/K - keep in mechanical embedding, remove otherwise
   * Copy QM-QM only if requested to do QM-QM LJ
   * Copy QM-MM only if electrostatic or polarizable embedding
   * Never copy MM-MM
   */
  // Decide what to copy and what to erase
  bool copy = false;
  bool erase = false;
  auto decide_copy_erase = [&](unsigned i, unsigned j)->void {
    copy = false;
    erase = false;
    switch (bool(topo.is_qm(i)) + bool(topo.is_qm(j))) {
      case 0 : { // MM-MM
        break;
      }
      case 1 : { // QM-MM
        // if MM is in QM buffer zone
        if (topo.is_qm_buffer(i) || topo.is_qm_buffer(j)) {
          copy = sim.param().qmmm.qm_lj;
          erase = true;
        }
        else if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
          copy = erase = true;
        }
        break;
      }
      case 2 : { // QM-QM
        copy = sim.param().qmmm.qm_lj;
        erase = true;
        break;
      }
    }
    return;
  };
  
  for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i) {
    const bool i_is_qm = topo.is_qm(i);
    if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
        && !i_is_qm && !use_qm_buffer) continue;

    // One-four pairs - they use CS6 and CS12
    for (topology::excl_cont_t::value_type::const_iterator
          it = topo.one_four_pair(i).begin()
        ; it != topo.one_four_pair(i).end(); ) {
      decide_copy_erase(i, *it);
      if (copy) {
        DEBUG(9, "Making copy of 1,4: " << i << " - " << *it);
        topo.qm_one_four_pair(i).insert(*it);
        topo.qm_all_exclusion(i).insert(*it);
      }
      if (erase) {
        DEBUG(9, "Removing 1,4 pair: " << i << " - " << *it);
        it = topo.one_four_pair(i).erase(it);
        continue;
      }
      ++it;
    }

    // Exclusions
    for (topology::excl_cont_t::value_type::const_iterator
          it = topo.exclusion(i).begin()
        ; it != topo.exclusion(i).end(); ) {
      decide_copy_erase(i, *it);
      if (copy) {
        DEBUG(9, "Making copy of exclusion: " << i << " - " << *it);
        topo.qm_exclusion(i).insert(*it);
        topo.qm_all_exclusion(i).insert(*it);
      }
      if (erase) {
        DEBUG(9, "Removing exclusion: " << i << " - " << *it);
        it = topo.exclusion(i).erase(it);
        continue;
      }
      ++it;
    }
  }

  // LJ Exceptions - they use LJEX
  for (std::vector<topology::lj_exception_struct>::iterator
        it = topo.lj_exceptions().begin()
      ; it != topo.lj_exceptions().end(); ) {
    const unsigned i = it->i;
    const unsigned j = it->j;
    assert(i < j);
    decide_copy_erase(i, j);
    if (copy) {
      DEBUG(9, "Making copy of LJ exception: " << i << " - " << j);
      topo.qm_lj_exceptions().push_back(*it);
      topo.qm_all_exclusion(i).insert(j);
    }
    if (erase) {
      DEBUG(9, "Removing LJ exception: " << i << " - " << j);
      it = topo.lj_exceptions().erase(it);
      continue;
    }
    ++it;
  }
  topo.update_all_exclusion();
}

void interaction::QMMM_Interaction::store_set_data(
        const topology::Topology& topo,
        configuration::Configuration& conf,
        const simulation::Simulation& sim
        ) {
  m_timer.start("set data summation");
  std::vector<QMMM_Nonbonded_Set *>::iterator
    it = m_qmmm_nonbonded_set.begin(),
    to = m_qmmm_nonbonded_set.end();

  // add the forces, energies, virial...
  for (; it != to; ++it) {
    DEBUG(7, "adding forces from set " << it - m_qmmm_nonbonded_set.begin());
    (*it)->update_configuration(topo, conf, sim);
  }
  m_timer.stop("set data summation");
}

int interaction::QMMM_Interaction::print_pairlist(const topology::Topology& topo
                                                , std::ostream & os) {
  DEBUG(4, "QMMM_Interaction::print_pairlist");

  Pairlist temp_solute, temp_solvent;
  temp_solute.resize(topo.num_atoms());
  temp_solvent.resize(topo.num_atoms());

  for (unsigned int atom_i = 0; atom_i < topo.num_atoms(); ++atom_i) {

    for (unsigned i = 0; i < m_set_size; ++i) {

      assert(m_qmmm_nonbonded_set.size() > unsigned(i));
      assert(m_qmmm_nonbonded_set[i]->pairlist().solute_short.size() > atom_i);
      assert(m_qmmm_nonbonded_set[i]->pairlist().solvent_short.size() > atom_i);

      for (unsigned int atom_j = 0;
              atom_j < m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i].size();
              ++atom_j) {

        assert(temp_solute.size() > atom_i);
        assert(temp_solute.size() > m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);

        if (m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j] < atom_i)
          temp_solute[m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]].push_back(atom_i);
        else
          temp_solute[atom_i].push_back(m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);
      }
      for (unsigned int atom_j = 0;
              atom_j < m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i].size();
              ++atom_j) {

        assert(temp_solvent.size() > atom_i);
        assert(temp_solvent.size() > m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);

        if (m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j] < atom_i)
          temp_solvent[m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]].push_back(atom_i);
        else
          temp_solvent[atom_i].push_back(m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);
      }
    }
  }
  os << temp_solute << std::endl
     << temp_solvent << std::endl;
  return 0;
}

void interaction::QMMM_Interaction::print_timing(std::ostream & os)
{
  if (m_rank == 0) {
    m_timer.print(os);
    m_worker->timer().print(os);
  }
}
