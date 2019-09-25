/**
 * @file local_elevation_interaction.cc
 * apply LEUS
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/special/qmmm/qm_storage.h"
#include "../../interaction/special/qmmm/mm_atom.h"
#include "../../interaction/special/qmmm/qm_worker.h"
#include "../../interaction/special/qmmm/mndo_worker.h"
#include "../../interaction/special/qmmm/mopac_worker.h"
#include "../../interaction/special/qmmm_interaction.h"
#include "../../interaction/special/qmmm/gathering.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"
#include "../../math/boundary_checks.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special
#define MAXPATH 10240

interaction::QMMM_Interaction::~QMMM_Interaction() {
  if (worker != NULL)
    delete worker;
}

configuration::Configuration interaction::QMMM_Interaction::AddRemove(topology::Topology &topo,
                                                                      configuration::Configuration &conf,
                                                                      simulation::Simulation &sim) {
    /*create a conf for capped system. Required for the evaluation of the forces on the external
    point charges in the MOPAC worker
    */
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

int interaction::QMMM_Interaction::prepare(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim) {
  // determine which MM atoms to include as point charges
  m_timer.start("MM atoms determination");
  interaction::determine_mm_atoms(topo, conf, sim, mm_atoms);
  m_timer.stop("MM atoms determination");

  // get the position of the QM atoms and gather
  m_timer.start("gathering");
  {
    qm_pos.resize(topo.qm_zone().size());
    unsigned int i = 0;
    for (std::set<topology::qm_atom_struct>::const_iterator
      it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++i) {
      qm_pos(i) = conf.current().pos(it->index);
    }
  }
  interaction::gather_qmzone(topo, conf, sim, qm_pos);
  interaction::gather_mm_atoms(topo, conf, sim, qm_pos, mm_atoms);
  m_timer.stop("gathering");
   
  return 0;
}

int interaction::QMMM_Interaction::
calculate_interactions(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim) {
  m_timer.start();
  storage.zero();
  
  // if polarisation is in use this was already called.
  if (!sim.param().polarise.cos) 
    prepare(topo, conf, sim);

  configuration::Configuration qmmm_conf=  conf;
  if ( topo.qm_mm_pair().size()  > 0 ) {
     configuration::Configuration qmmm_conf=  interaction::QMMM_Interaction::AddRemove(topo,conf,sim);
  }
/*
    std::cout << "conf corrent at step "  << sim.steps() << std::endl;
       for (unsigned int i = 0; i < 26; i++) {
           std::cout << "conf corrent force before:   " << i << "  " << conf.current().force(i)[0] <<
                     "   " << conf.current().force(i)[1] << "  "
                     << conf.current().force(i)[2] << std::endl;
       }
       */
  if (worker->run_QM(topo, conf, sim, qm_pos, mm_atoms, storage,stor_link,qmmm_conf)) {
    m_timer.stop();
    return 1;
  }

  // add QM forces
  AddRemove2(topo,conf,sim);
  for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
    //std::cout<< i << "  storage.force: " << storage.force(i)[0] << " " <<
   //          storage.force(i)[1] << "  " << storage.force(i)[2] << std::endl;
    conf.current().force(i) += storage.force(i);
    if (topo.is_polarisable(i)) {
      conf.current().force(i) += storage.cos_force(i);
    }
  }
  
  conf.current().energies.qm_total = storage.energy;
  m_timer.stop();
  return 0;
}

int interaction::QMMM_Interaction::add_electric_field_contribution(topology::Topology& topo,
        configuration::Configuration& conf, simulation::Simulation& sim, 
        math::VArray & electric_field) {
  
  m_timer.start();
  m_timer.start("polarisation");
  
  storage.zero();

  if (worker->run_QM(topo, conf, sim, qm_pos, mm_atoms, storage,stor_link,conf)) {
    m_timer.stop();
    m_timer.stop("polarisation");
    return 1;
  }
  for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
    if (!topo.is_polarisable(i))
      continue;

    // get electric field at either the charge or the COS site.
    math::Vec e;
    switch (sim.param().polarise.efield_site) {
      case simulation::ef_atom:
        e = storage.force(i) / (topo.charge(i) - topo.coscharge(i));
        break;
      case simulation::ef_cos:
        e = storage.cos_force(i) / topo.coscharge(i);
        break;
      default:
        io::messages.add("Electric field site not implemented.",
                "QMMM_Interaction", io::message::critical);
        m_timer.stop("polarisation");
        m_timer.stop();
        return 1;
    }
    // add the electric field contribution
    electric_field(i) += e;
  }
  
  m_timer.stop("polarisation");
  m_timer.stop();
  return 0;
}
int interaction::QM_Worker::qm_idGenerator=1;

int interaction::QMMM_Interaction::init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os,
            bool quiet) { 

    quiet = false;
    qmmm_topo = new topology::Topology(topo, 1); //copy topology in order to restore MM-terms
    //required for AddRemove scheme
  if (topo.qm_zone().empty()) {
    io::messages.add("No QM zone defined", "QMMM_Interaction", io::message::error);
    return 1;
  }
    stor_link.resize(topo.qm_mm_pair().size(),topo.qm_mm_pair().size());
  storage.resize(topo.num_atoms(), topo.qm_zone().size());
  
  worker = interaction::QM_Worker::get_instance(sim,1);
  if (worker == NULL) return 1;
  if (worker->init(topo, conf, sim)) return 1;
  
  if (!quiet) {
    os << "QM/MM" << std::endl
            << "  QM worker: " << worker->name() << std::endl;

    os << "  QM zone: " << std::endl;
    for (std::set<topology::qm_atom_struct>::const_iterator
      it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
            os << "     - " << std::setw(5) << it->index + 1 << std::setw(5) << it->atomic_number;
      if (it->link) {
        os << " link atom";
      }
      os << std::endl;
    }
    os.precision(3);
    if (sim.param().qmmm.cutoff > 0.0) {
      os << "  Cutoff for inclusion of MM atoms: " << std::setw(8)
              << sim.param().qmmm.cutoff << std::endl;
    } else {
      os << "  All atoms are included as MM atoms." << std::endl;
    }
    if (sim.param().qmmm.mmscal > 0) {
      os << "  MM point charges are scaled with 2/pi atan(s*|R|) with s=  " <<
         sim.param().qmmm.mmscal << std::endl;
    }
    os << "END" << std::endl;
  }
  
  if (sim.param().qmmm.cutoff > 0.0 && 
          !math::boundary_check_cutoff(conf.current().box, conf.boundary_type,
          sim.param().qmmm.cutoff)) {
    io::messages.add("The cutoff RCUTQ is too large for the size of the "
            "computational box.", "QMMM_Interaction", io::message::error);
    return 1;
  }
  
  if (sim.param().qmmm.cutoff == 0.0 && conf.boundary_type != math::vacuum) {
    io::messages.add("A cutoff RCUTQ >= 0.0 has to be used for simulations "
            "using periodic boundary conditions.", "QMMM_Interaction", 
            io::message::error);
    return 1;
  }
  
  // set charge of QM atoms to 0
  for(std::set<topology::qm_atom_struct>::const_iterator
      it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
    topo.charge()[it->index] = 0.0;
  }
    //also set charges of link atoms to zero
    //note that bonds should only been cut between two neutral atoms anyway...
    for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
                 it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it){
    //    std::cout << "it->second link atoms to zero " << it->second << std::endl;
        topo.charge()[it->first] = 0.0;
        topo.charge()[it->second] = 0.0;
}
    if (false) {
        for (unsigned int i = 0; i < topo.solute().bonds().size(); ++i) {
            std::cout << "i  bonds " << topo.solute().bonds()[i].i << std::endl;
            for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
                         it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it) {
                std::cout << "it->second " << it->second << " bonds " << topo.solute().bonds()[i].i <<
                          topo.solute().bonds()[i].j << "  " << std::endl;
                if (topo.solute().bonds()[i].i == it->second) {
                    std::cout << "setting charge of " << topo.solute().bonds()[i].j << "  to 0.0" << std::endl;
                } else if (topo.solute().bonds()[i].j == it->second) {
                    std::cout << "setting charge of " << topo.solute().bonds()[i].i << "  to 0.0" << std::endl;
                }
            }
        }
    }
    /*
  topology::Chargegroup_Iterator cg = topo.chargegroup_begin();
  topology::Chargegroup_Iterator to = topo.chargegroup_end();
    for(unsigned int cg = 0; cg < topo.num_solute_chargegroups(); ++cg) {
      DEBUG(10, "Initial QM-charges: " <<     cg << " charge: " << topo.charge()[cg] );
  
    }
     */
  // remove bonded terms
    //now we only remove bonded contributions where all atoms are in the QM-zone
    //se Walker et. al J Comput Chem 29 2007
  for (unsigned int i = 0; i < topo.solute().bonds().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().bonds()[i].i) &&
        topo.in_qm_zone(topo.solute().bonds()[i].j)) {
      topo.solute().bonds().erase(topo.solute().bonds().begin() + i);
  std::cout << "i:" << i  << "erased bonds: " <<
          topo.solute().bonds()[i].i << " " << topo.solute().bonds()[i].j << std::endl;
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().cgbonds().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().cgbonds()[i].i) || 
        topo.in_qm_zone(topo.solute().cgbonds()[i].j)) {
      topo.solute().cgbonds().erase(topo.solute().cgbonds().begin() + i);
       std::cout << "i:" << i  << "erased cgbonds: " << 
              topo.solute().cgbonds()[i].j << std::endl;
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().distance_constraints().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().distance_constraints()[i].i) &&
        topo.in_qm_zone(topo.solute().distance_constraints()[i].j)) {
      topo.solute().distance_constraints().erase(topo.solute().distance_constraints().begin() + i);
        /*
      std::cout << "i:" << i  << "erased distance constraint: " <<
              topo.solute().distance_constraints()[i].i << std::endl;
      std::cout << "j:" << i  << "erased distance constraint: " << 
              topo.solute().distance_constraints()[i].j << std::endl;
      std::cout << "type:" << i  << "erased distance constraint: " << 
              topo.solute().distance_constraints()[i].type << std::endl;
        */
      --i;
    }
  }
  
  
  for (unsigned int i = 0; i < topo.solute().angles().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().angles()[i].i) &&
        topo.in_qm_zone(topo.solute().angles()[i].j) &&
        topo.in_qm_zone(topo.solute().angles()[i].k)) {
        std::cout << "i:" << i  << "erased angles: " <<
              topo.solute().angles()[i].i << "     " <<
              topo.solute().angles()[i].j << "     " << 
              topo.solute().angles()[i].k << std::endl;
      topo.solute().angles().erase(topo.solute().angles().begin() + i);
      --i;
    }
  }
    std::vector<topology::four_body_term_struct> lost_diheds;
  for (unsigned int i = 0; i < topo.solute().dihedrals().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().dihedrals()[i].i) || 
        topo.in_qm_zone(topo.solute().dihedrals()[i].j) || 
        topo.in_qm_zone(topo.solute().dihedrals()[i].k) ||
        topo.in_qm_zone(topo.solute().dihedrals()[i].l)) {
        lost_diheds.push_back(topo.solute().dihedrals()[i]);
      topo.solute().dihedrals().erase(topo.solute().dihedrals().begin() + i);
      --i;
    }
  }
    unsigned  int i =0;
    for (std::vector<topology::four_body_term_struct>::const_iterator
                 it2 = lost_diheds.begin(),
                 to2 = lost_diheds.end(); it2 != to2; ++it2, i++) {
        std::cout << "i:" << i << "erased dihedrals: "
                  << it2->i + 1 << "   "
                  << it2->j + 1 << "   "
                  << it2->k + 1 << "   "
                  << it2->l + 1 << "   " << std::endl;
    }
  for (unsigned int i = 0; i < topo.solute().improper_dihedrals().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().improper_dihedrals()[i].i) &&
        topo.in_qm_zone(topo.solute().improper_dihedrals()[i].j) &&
        topo.in_qm_zone(topo.solute().improper_dihedrals()[i].k) &&
        topo.in_qm_zone(topo.solute().improper_dihedrals()[i].l)) {
        /*
        std::cout << "i:" << i  << "erased improper dihedrals: " <<
              topo.solute().improper_dihedrals()[i].i << "    " <<
              topo.solute().improper_dihedrals()[i].j << "    " <<
              topo.solute().improper_dihedrals()[i].k << "    " <<
              topo.solute().improper_dihedrals()[i].l << std::endl;
              */
      topo.solute().improper_dihedrals().erase(topo.solute().improper_dihedrals().begin() + i);
      --i;
    }
  }
  for (unsigned int i = 0; i < topo.solute().crossdihedrals().size(); ++i) {
    if (topo.in_qm_zone(topo.solute().crossdihedrals()[i].a) &&
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].b) &&
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].c) &&
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].d) &&
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].e) &&
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].f) &&
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].g) &&
        topo.in_qm_zone(topo.solute().crossdihedrals()[i].h)) {
        std::cout << "i:" << i  << "erased crossdihedrals: " <<
              topo.solute().crossdihedrals()[i].a << std::endl;
      topo.solute().crossdihedrals().erase(topo.solute().crossdihedrals().begin() + i);
      --i;
    }
  }
  if (false ) {
      //restoring the bonded terms, where the MM-Link atom is involved
      unsigned int m1, q1;
      for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
                   it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it) {
          q1 = it->first;
          m1 = it->second;
          unsigned int i = 0;

          //bonds
          for (std::vector<topology::two_body_term_struct>::const_iterator
                       it2 = qmmm_topo->solute().bonds().begin(),
                       to2 = qmmm_topo->solute().bonds().end(); it2 != to2; ++it2, i++) {
              if (it2->i == m1 || it2->j == m1) {
                  if (topo.in_qm_zone(it2->i) == 0 && topo.in_qm_zone(it2->j) == 0) continue;

                  else {
                      topo.solute().bonds().push_back(
                              topology::two_body_term_struct(qmmm_topo->solute().bonds()[i].i,
                                                             qmmm_topo->solute().bonds()[i].j,
                                                             qmmm_topo->solute().bonds()[i].type));
                      std::cout << "i:" << i << "QMMM TOPO  restored bonds: " <<
                                qmmm_topo->solute().bonds()[i].i + 1 << "   " <<
                                qmmm_topo->solute().bonds()[i].j + 1 << std::endl;
//                 sim.param().qmmm.QM_CAP_topo->solute().push
                  }
              }
          }

          //angles
          i = 0;
          for (std::vector<topology::three_body_term_struct>::const_iterator
                       it2 = qmmm_topo->solute().angles().begin(),
                       to2 = qmmm_topo->solute().angles().end(); it2 != to2; ++it2, i++) {
              if (it2->i == m1 || it2->j == m1 || it2->k == m1) {
                  unsigned int mm_count = 0;
                  if (topo.in_qm_zone(it2->i) == 0) mm_count++;
                  if (topo.in_qm_zone(it2->j) == 0) mm_count++;
                  if (topo.in_qm_zone(it2->k) == 0) mm_count++;

                  if (mm_count > 1) continue;
                  else {
                      topo.solute().angles().push_back(
                              topology::three_body_term_struct(qmmm_topo->solute().angles()[i].i,
                                                               qmmm_topo->solute().angles()[i].j,
                                                               qmmm_topo->solute().angles()[i].k,
                                                               qmmm_topo->solute().angles()[i].type));

                      std::cout << "i:" << i << "QMMM TOPO  restored angles: " <<
                                qmmm_topo->solute().angles()[i].i + 1 << "     " <<
                                qmmm_topo->solute().angles()[i].j + 1 << "     " <<
                                qmmm_topo->solute().angles()[i].k + 1 << ::std::endl;
                  }
              }
          }
          //distance constraints
          i = 0;
          for (std::vector<topology::two_body_term_struct>::const_iterator
                       it2 = qmmm_topo->solute().distance_constraints().begin(),
                       to2 = qmmm_topo->solute().distance_constraints().end();
               it2 != to2; ++it2, i++) {
              if (it2->i == m1 || it2->j == m1) {
                  if ((topo.in_qm_zone(it2->i)) == 0 && topo.in_qm_zone(it2->j) == 0) continue;
                      //  if (topo.in_qm_zone(it2->j) == 0 && it2->j != m1) continue;

                  else {
                      topo.solute().distance_constraints().push_back(
                              topology::two_body_term_struct(
                                      qmmm_topo->solute().distance_constraints()[i].i,
                                      qmmm_topo->solute().distance_constraints()[i].j,
                                      qmmm_topo->solute().distance_constraints()[i].type));
                      std::cout << "i:" << i << "QMMM TOPO  restored distance constraints: " <<
                                //                 std::cout << "i:" << i  << "QMMM TOPO  removing distrance constraints: " <<
                                qmmm_topo->solute().distance_constraints()[i].i << "   " <<
                                qmmm_topo->solute().distance_constraints()[i].j << std::endl;

                  }
              }
          }
      }
      }
      /*
       dihedrals  [i,j,k,l]
       we want to describe the torsion around the MM-QM Atom classically
       consider the case [qm,qm,mm,mm] or [mm,mm,qm,qm]
      */
      i = 0;
      for (std::vector<topology::four_body_term_struct>::const_iterator
                   it2 = qmmm_topo->solute().dihedrals().begin(),
                   to2 = qmmm_topo->solute().dihedrals().end(); it2 != to2; ++it2, i++) {
          if (((topo.in_qm_zone(it2->i) + topo.in_qm_zone(it2->j) == 2) &&
               (topo.in_qm_zone(it2->k) + topo.in_qm_zone(it2->l) == 0)) || (
                      (topo.in_qm_zone(it2->i) + topo.in_qm_zone(it2->j) == 0) &&
                      (topo.in_qm_zone(it2->k) + topo.in_qm_zone(it2->l) == 2))
                  ) {
              topo.solute().dihedrals().push_back(
                      topology::four_body_term_struct(qmmm_topo->solute().dihedrals()[i].i,
                                                      qmmm_topo->solute().dihedrals()[i].j,
                                                      qmmm_topo->solute().dihedrals()[i].k,
                                                      qmmm_topo->solute().dihedrals()[i].l,
                                                      qmmm_topo->solute().dihedrals()[i].type));
              std::cout << "i:" << i << "QMMM TOPO  restored dihedrals: " <<
                        qmmm_topo->solute().dihedrals()[i].i + 1 << "     " <<
                        qmmm_topo->solute().dihedrals()[i].j + 1 << "     " <<
                        qmmm_topo->solute().dihedrals()[i].k + 1 << "     " <<
                        qmmm_topo->solute().dihedrals()[i].l + 1 << ::std::endl;
          }
      }

      /*
        dihedrals  [i,j,k,l]
      we want to describe the torsion around the MM-QM Atom classically
      consider the case [qm,mm,mm,mm]
      */
      i = 0;
      for (std::vector<topology::four_body_term_struct>::const_iterator
                   it2 = qmmm_topo->solute().dihedrals().begin(),
                   to2 = qmmm_topo->solute().dihedrals().end(); it2 != to2; ++it2, i++) {
          if (topo.in_qm_zone(it2->i) == 1 && topo.in_qm_zone(it2->l) == 0)
              if (topo.in_qm_zone(it2->k) + topo.in_qm_zone(it2->j) == 0) {
                  topo.solute().dihedrals().push_back(
                          topology::four_body_term_struct(qmmm_topo->solute().dihedrals()[i].i,
                                                          qmmm_topo->solute().dihedrals()[i].j,
                                                          qmmm_topo->solute().dihedrals()[i].k,
                                                          qmmm_topo->solute().dihedrals()[i].l,
                                                          qmmm_topo->solute().dihedrals()[i].type));
                  std::cout << "i:" << i << "QMMM TOPO  restored dihedrals: " <<
                            qmmm_topo->solute().dihedrals()[i].i + 1 << "     " <<
                            qmmm_topo->solute().dihedrals()[i].j + 1 << "     " <<
                            qmmm_topo->solute().dihedrals()[i].k + 1 << "     " <<
                            qmmm_topo->solute().dihedrals()[i].l + 1 << ::std::endl;
              }
      }
      /*
      dihedrals  [i,j,k,l]
      we want to describe the torsion around the MM-QM Atom classicaly
      consider the case [mm,mm,mm,qm]
      */
      i = 0;
      for (std::vector<topology::four_body_term_struct>::const_iterator
                   it2 = qmmm_topo->solute().dihedrals().begin(),
                   to2 = qmmm_topo->solute().dihedrals().end(); it2 != to2; ++it2, i++) {
          if (topo.in_qm_zone(it2->i) == 0 && topo.in_qm_zone(it2->l) == 1)
              if (topo.in_qm_zone(it2->k) + topo.in_qm_zone(it2->j) == 0) {
                  topo.solute().dihedrals().push_back(
                          topology::four_body_term_struct(qmmm_topo->solute().dihedrals()[i].i,
                                                          qmmm_topo->solute().dihedrals()[i].j,
                                                          qmmm_topo->solute().dihedrals()[i].k,
                                                          qmmm_topo->solute().dihedrals()[i].l,
                                                          qmmm_topo->solute().dihedrals()[i].type));
                  std::cout << "i:" << i << "QMMM TOPO  restored dihedrals: " <<
                            qmmm_topo->solute().dihedrals()[i].i + 1 << "     " <<
                            qmmm_topo->solute().dihedrals()[i].j + 1 << "     " <<
                            qmmm_topo->solute().dihedrals()[i].k + 1 << "     " <<
                            qmmm_topo->solute().dihedrals()[i].l + 1 << ::std::endl;
              }
      }

/*
      //improper dihedrals
      //fixme has to be updated, but at the moment we should not cut through such bonds anyway.
      i = 0;
      for (std::vector<topology::four_body_term_struct>::const_iterator
                   it2 = qmmm_topo->solute().improper_dihedrals().begin(),
                   to2 = qmmm_topo->solute().improper_dihedrals().end(); it2 != to2; ++it2, i++) {
          if (it2->i == m1 || it2->j == m1 || it2->k == m1 || it2->l == m1) {
              unsigned int mm_count = 0;
              if (topo.in_qm_zone(it2->i) == 0) mm_count++;
              if (topo.in_qm_zone(it2->j) == 0) mm_count++;
              if (topo.in_qm_zone(it2->k) == 0) mm_count++;
              if (topo.in_qm_zone(it2->l) == 0) mm_count++;

              if (mm_count > 3) continue;
              else {
                  topo.solute().improper_dihedrals().push_back(
                          topology::four_body_term_struct(
                                  qmmm_topo->solute().improper_dihedrals()[i].i,
                                  qmmm_topo->solute().improper_dihedrals()[i].j,
                                  qmmm_topo->solute().improper_dihedrals()[i].k,
                                  qmmm_topo->solute().improper_dihedrals()[i].l,
                                  qmmm_topo->solute().improper_dihedrals()[i].type));

                  std::cout << "i:" << i << "QMMM TOPO  restored improper_dihedrals: " <<
                            qmmm_topo->solute().improper_dihedrals()[i].i + 1 << "     " <<
                            qmmm_topo->solute().improper_dihedrals()[i].j + 1 << "     " <<
                            qmmm_topo->solute().improper_dihedrals()[i].k + 1 << "     " <<
                            qmmm_topo->solute().improper_dihedrals()[i].l + 1 << ::std::endl;
              }
          }
      }
      */
     // Also exclude non-bonded within QM zone if asked to

  return 0;
}

