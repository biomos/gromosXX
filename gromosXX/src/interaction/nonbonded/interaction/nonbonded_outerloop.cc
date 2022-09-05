/**
 * @file nonbonded_outerloop.cc
 * template methods of Nonbonded_Outerloop.
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
#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/latticesum.h"
#include "../../../configuration/mesh.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_innerloop.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"

#include "../../../interaction/nonbonded/interaction_spec.h"

#include "../../../util/debug.h"
#include "../../../interaction/nonbonded/innerloop_template.h"

#include "../../../simulation/parameter.h"

#include "../../../interaction/qmmm/qmmm_interaction.h"

#include "../../interaction.h"

#ifdef OMP
#include <omp.h>
#endif

#include <algorithm>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

#ifdef OMP
    math::VArray interaction::Nonbonded_Outerloop::electric_field = 0.0;
    double interaction::Nonbonded_Outerloop::minfield = 0.0;
#endif

/**
 * Constructor.
 */
interaction::Nonbonded_Outerloop
::Nonbonded_Outerloop(Nonbonded_Parameter &nbp)
: m_param(nbp) {
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
        Storage & storage,
        bool longrange, util::Algorithm_Timer & timer, bool master) {
  SPLIT_INNERLOOP(_lj_crf_outerloop, topo, conf, sim,
          pairlist_solute, pairlist_solvent, storage, longrange, timer, master);
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
        Storage & storage, bool longrange,
        util::Algorithm_Timer & timer, bool master) {
  DEBUG(7, "\tcalculate interactions");

  // WORKAROUND! See definition of _lj_crf_outerloop_fast
  if (t_interaction_spec::boundary_type == math::rectangular &&
      t_interaction_spec::interaction_func == simulation::lj_crf_func &&
      sim.param().innerloop.method != simulation::sla_cuda) {
    _lj_crf_outerloop_fast<t_interaction_spec::charge_type>(topo, conf, sim,
                        pairlist_solute, pairlist_solvent,
                        storage, longrange, timer, master);
    return;
  }

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  /*
    variables for a OMP parallelizable loop.
    outer index has to be integer...
   */
  std::vector<unsigned int>::const_iterator j_it, j_to;

  unsigned int size_i = unsigned(pairlist_solute.size());
  DEBUG(10, "lj_crf_outerloop pairlist size " << size_i);

  const unsigned int end = topo.num_solute_atoms();

  unsigned int i = 0;
  for (i = 0; i < end; ++i) {
    for (j_it = pairlist_solute[i].begin(),
            j_to = pairlist_solute[i].end();
            j_it != j_to;
            ++j_it) {

      DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << *j_it);

      // shortrange, therefore store in simulation.system()
      innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }
/*only for DEBUG*/
 /* DEBUG(1,"current solute pairlist:\n");
unsigned int i_deb;
  for (i_deb=0; i_deb < end; ++i_deb) {
    for (j_it = pairlist_solute[i_deb].begin(),
            j_to = pairlist_solute[i_deb].end();
            j_it < j_to;
            j_it++) {
      //DEBUG(10, "i " << i_deb << " j " << *j_it << " i " << topo.solvent(0).atom(i_deb).name << " j " << topo.solvent(0).atom(*j_it).name);
      DEBUG(1, "i " << i_deb << " j " << *j_it);
    }
  }
  DEBUG(1,"current solvent pairlist:\n");
  for (; i_deb < size_i; ++i_deb) {
    for (j_it = pairlist_solvent[i_deb].begin(),
            j_to = pairlist_solvent[i_deb].end();
            j_it < j_to;
            j_it++) {
      //DEBUG(10, "i " << i_deb << " j " << *j_it << " i " << topo.solvent(0).atom(i_deb).name << " j " << topo.solvent(0).atom(*j_it).name);
      DEBUG(1, "i " << i_deb << " j " << *j_it);
    }
  }
  */
/*DEBUG end*/

  // cuda doesn't do solvent-solvent here
  if (sim.param().innerloop.method == simulation::sla_cuda) return;
  // solvent-solvent
  const std::string timer_name(longrange ? "longrange solvent-solvent" : "solvent-solvent");
  if (master)
    timer.start(timer_name);
  if (sim.param().force.special_loop == simulation::special_loop_spc) { // special solvent loop
    // solvent - solvent with spc innerloop...
    for (; i < size_i; i += 3) { // use every third pairlist (OW's)
      for (j_it = pairlist_solvent[i].begin(),
              j_to = pairlist_solvent[i].end();
              j_it != j_to;
              j_it += 3) { // use every third atom (OW) in pairlist i

        DEBUG(10, "\tsolvent-solvent longrange spc_nonbonded_interaction: i " << i << " j " << *j_it);

        innerloop.spc_innerloop(topo, conf, i, *j_it, storage, periodicity);
      }
    } 
  } else if (sim.param().force.special_loop == simulation::special_loop_spc_table) { // special solvent loop
    // solvent - solvent with tabulated spc innerloop...
    if (longrange) {
      // here we call the longrange function that uses a smaller table
      for (; i < size_i; i += 3) { // use every third pairlist (OW's)
        for (j_it = pairlist_solvent[i].begin(),
                j_to = pairlist_solvent[i].end();
                j_it != j_to;
                j_it += 3) { // use every third atom (OW) in pairlist i

          DEBUG(10, "\tsolvent-solvent (tabulated) longrange spc_nonbonded_interaction: i " << i << " j " << *j_it);
          //innerloop.spc_innerloop(topo, conf, i, *j_it, storage, periodicity);
          innerloop.longrange_spc_table_innerloop(topo, conf, i, *j_it, storage, periodicity);
        }
      }
    } else { // shortrange
      for (; i < size_i; i += 3) { // use every third pairlist (OW's)
        for (j_it = pairlist_solvent[i].begin(),
                j_to = pairlist_solvent[i].end();
                j_it != j_to;
                j_it += 3) { // use every third atom (OW) in pairlist i

          DEBUG(10, "\tsolvent-solvent (tabulated) shortrange spc_nonbonded_interaction: i " << i << " j " << *j_it);
          //innerloop.spc_innerloop(topo, conf, i, *j_it, storage, periodicity);
          innerloop.shortrange_spc_table_innerloop(topo, conf, i, *j_it, storage, periodicity);
        }
      }
    }
  } else if (sim.param().force.special_loop == simulation::special_loop_generic) {
    const unsigned int num_solvent_atoms = topo.solvent(0).num_atoms();
    // prepare parameters
    const unsigned int num_param = num_solvent_atoms * num_solvent_atoms;
    typename interaction::Nonbonded_Innerloop<t_interaction_spec>::solvent_pair_parameter pair_parameter[num_param];

    unsigned int param = 0;
    for (unsigned int atom_i = 0; atom_i < num_solvent_atoms; ++atom_i) {
      for (unsigned int atom_j = 0; atom_j < num_solvent_atoms; ++atom_j, ++param) {
        assert(param < num_param);
        DEBUG(10, "\tsolvent pair parameter: " << param);

        const lj_parameter_struct & lj =
                innerloop.param()->lj_parameter(topo.solvent(0).atom(atom_i).iac, topo.solvent(0).atom(atom_j).iac);
        pair_parameter[param].c12 = lj.c12;
        pair_parameter[param].c6 = lj.c6;

        pair_parameter[param].q = math::four_pi_eps_i *
                topo.solvent(0).atom(atom_i).charge *
                topo.solvent(0).atom(atom_j).charge;

        DEBUG(10, "\t\tc12: " << pair_parameter[param].c12 << " c6: " <<
                pair_parameter[param].c6 << " q: " << pair_parameter[param].q);
      }
    }

    // use num_solvent_atoms-th atom (first of solvent molecule i)
    for (; i < size_i; i += num_solvent_atoms) {
      for (j_it = pairlist_solvent[i].begin(),
              j_to = pairlist_solvent[i].end();
              j_it != j_to;
              j_it += num_solvent_atoms) { // use num_solvent_atoms-th atom (first of solvent molecule j)

        DEBUG(10, "\tsolvent_nonbonded_interaction (special_loop_generic): i " << i << " j " << *j_it);

        innerloop.solvent_innerloop(topo, pair_parameter, conf, num_solvent_atoms, i, *j_it, storage, periodicity);
      }
    }
  } else { // normal solvent loop
    for (; i < size_i; ++i) {
      for (j_it = pairlist_solvent[i].begin(),
              j_to = pairlist_solvent[i].end();
              j_it != j_to;
              ++j_it) {

        DEBUG(10, "\tsolvent_nonbonded_interaction (normal solvent loop): i " << i << " j " << *j_it);

        innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
      }
    }
  }
  if (master)
    timer.stop(timer_name);
}


/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
// WORKAROUND - see definition!
template<simulation::charge_type_enum t_charge_type>
void interaction::Nonbonded_Outerloop
::_lj_crf_outerloop_fast(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Pairlist const & pairlist_solute,
        Pairlist const & pairlist_solvent,
        Storage & storage, bool longrange,
        util::Algorithm_Timer & timer, bool master) {
  DEBUG(7, "\tcalculate interactions lj_crf_outerloop_fast");

  math::Periodicity<math::rectangular> periodicity(conf.current().box);
  periodicity.recalc_shift_vectors();
  Nonbonded_Innerloop<interaction::Interaction_Spec<math::rectangular,
                      simulation::lj_crf_func,
                      t_charge_type> > innerloop(m_param);
  innerloop.init(sim);

  /*
    variables for a OMP parallelizable loop.
    outer index has to be integer...
   */
  std::vector<unsigned int>::const_iterator j_it, j_to;

  unsigned int size_i = unsigned(pairlist_solute.size());
  DEBUG(10, "lj_crf2 outerloop pairlist size " << size_i);

  const unsigned int end = topo.num_solute_atoms();

  unsigned int i = 0;
  for (i = 0; i < end; ++i) {
    const math::Vec posI = conf.current().pos(i);
    const unsigned int eg_i = topo.atom_energy_group(i);
    math::Vec groupForce(0.0);
    int k = 0;
    

    // shortrange, therefore store in simulation.system()
    for (j_it = pairlist_solute[i].begin(),
         j_to = pairlist_solute[i].end();
         j_it != j_to;
         ++j_it) {
      DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << *j_it);
            
      math::Vec r;
      const int kk = periodicity.nearest_image(posI, conf.current().pos(*j_it), r);
      
      if (kk != k) {
        storage.force(i) += groupForce;
        if (k != 0) {
          const math::Vec shift = periodicity.shift(k+13).pos;
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
              storage.virial_tensor(b, a) += shift(b) * groupForce(a);
            }
          }
        }
        groupForce = 0.;
        k = kk;
        
      }
      
      const double dist2 = abs2(r);
      math::Vec force;
      double f = 0.0;
      double e_lj = 0.0, e_crf = 0.0;
      
      innerloop.lj_crf_innerloop_2(topo, i, *j_it, dist2, f, e_lj, e_crf);      
      
      const unsigned int eg_j = topo.atom_energy_group(*j_it);
      DEBUG(11, "\tenergy group i " << eg_i << " j " << eg_j);
      storage.energies.lj_energy[eg_i][eg_j] += e_lj;
      storage.energies.crf_energy[eg_i][eg_j] += e_crf;
            
      force = f * r;
      storage.force(*j_it) -= force;
      groupForce += force;
      
      // shortrange, therefore store in simulation.system()
      // innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
    
    storage.force(i) += groupForce;
    if (k != 0) {
      const math::Vec shift = periodicity.shift(k+13).pos;
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
          storage.virial_tensor(b, a) += shift(b) * groupForce(a);
        }
      }
    }
  }
    
  for (unsigned int ii = 0; ii < topo.num_solute_atoms(); ++ii) {
    const math::Vec pos = conf.current().pos(ii);
    const math::Vec force = storage.force(ii);
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
        storage.virial_tensor(b, a) += pos(b) * force(a);
      }
    }
  }
/*only for DEBUG*/
 /* DEBUG(1,"current solute pairlist:\n");
unsigned int i_deb;
  for (i_deb=0; i_deb < end; ++i_deb) {
    for (j_it = pairlist_solute[i_deb].begin(),
            j_to = pairlist_solute[i_deb].end();
            j_it < j_to;
            j_it++) {
      //DEBUG(10, "i " << i_deb << " j " << *j_it << " i " << topo.solvent(0).atom(i_deb).name << " j " << topo.solvent(0).atom(*j_it).name);
      DEBUG(1, "i " << i_deb << " j " << *j_it);
    }
  }
  DEBUG(1,"current solvent pairlist:\n");
  for (; i_deb < size_i; ++i_deb) {
    for (j_it = pairlist_solvent[i_deb].begin(),
            j_to = pairlist_solvent[i_deb].end();
            j_it < j_to;
            j_it++) {
      //DEBUG(10, "i " << i_deb << " j " << *j_it << " i " << topo.solvent(0).atom(i_deb).name << " j " << topo.solvent(0).atom(*j_it).name);
      DEBUG(1, "i " << i_deb << " j " << *j_it);
    }
  }
  */
/*DEBUG end*/

  // cuda doesn't do solvent-solvent here
  if (sim.param().innerloop.method == simulation::sla_cuda) return;
  // solvent-solvent
  const std::string timer_name(longrange ? "longrange solvent-solvent" : "solvent-solvent");
  if (master)
    timer.start(timer_name);
  if (sim.param().force.special_loop == simulation::special_loop_spc) { // special solvent loop
    // solvent - solvent with spc innerloop...

    // only one energy group
    const int egroup = topo.atom_energy_group(topo.num_solute_atoms());

    for (; i < size_i; i += 3) { // use every third pairlist (OW's)

      const math::Vec posI  = conf.current().pos(i);
      const math::Vec posI1 = conf.current().pos(i + 1);
      const math::Vec posI2 = conf.current().pos(i + 2);
      math::Vec groupForce0(0.0);
      math::Vec groupForce1(0.0);
      math::Vec groupForce2(0.0);
      int k = 0;
      math::Vec shift = periodicity.shift(k + 13).pos;
      double tx = shift(0), ty = shift(1), tz = shift(2);

      double dist6i = 0.0;
      double e_lj = 0., e_crf = 0.;
      double r2[9], r2i[9], ri[9], x[9], y[9], z[9], f[9], fx[9], fy[9], fz[9];
      math::Vec r;

      for (j_it = pairlist_solvent[i].begin(),
          j_to = pairlist_solvent[i].end();
          j_it != j_to;
          j_it += 3) { // use every third atom (OW) in pairlist i

        DEBUG(10, "\tsolvent-solvent spc_nonbonded_interaction: i " << i << " j " << *j_it);

        const int kk = periodicity.nearest_image(posI, conf.current().pos(*j_it), r);

        if (kk != k) {
          storage.force(i) += groupForce0;
          storage.force(i + 1) += groupForce1;
          storage.force(i + 2) += groupForce2;
          
          if (k != 0) {
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) {
                //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
                storage.virial_tensor(b, a) += shift(b) * (groupForce0(a) + groupForce1(a) + groupForce2(a));
              }
            }
          }
          
          groupForce0 = 0.;
          groupForce1 = 0.;
          groupForce2 = 0.;
          k = kk;
          shift = periodicity.shift(k + 13).pos;
          tx = shift(0);
          ty = shift(1);
          tz = shift(2);
        }

        math::Vec const * const pos_j = &conf.current().pos(*j_it);
        math::Vec * const force_j = &storage.force(*j_it);

        x[0] = r(0);
        y[0] = r(1);
        z[0] = r(2);

        x[1] = posI(0) - (*(pos_j + 1))(0) + tx;
        y[1] = posI(1) - (*(pos_j + 1))(1) + ty;
        z[1] = posI(2) - (*(pos_j + 1))(2) + tz;

        x[2] = posI(0) - (*(pos_j + 2))(0) + tx;
        y[2] = posI(1) - (*(pos_j + 2))(1) + ty;
        z[2] = posI(2) - (*(pos_j + 2))(2) + tz;

        x[3] = posI1(0) - (*(pos_j))(0) + tx;
        y[3] = posI1(1) - (*(pos_j))(1) + ty;
        z[3] = posI1(2) - (*(pos_j))(2) + tz;

        x[4] = posI2(0) - (*(pos_j))(0) + tx;
        y[4] = posI2(1) - (*(pos_j))(1) + ty;
        z[4] = posI2(2) - (*(pos_j))(2) + tz;

        x[5] = posI1(0) - (*(pos_j + 1))(0) + tx;
        y[5] = posI1(1) - (*(pos_j + 1))(1) + ty;
        z[5] = posI1(2) - (*(pos_j + 1))(2) + tz;

        x[6] = posI1(0) - (*(pos_j + 2))(0) + tx;
        y[6] = posI1(1) - (*(pos_j + 2))(1) + ty;
        z[6] = posI1(2) - (*(pos_j + 2))(2) + tz;

        x[7] = posI2(0) - (*(pos_j + 1))(0) + tx;
        y[7] = posI2(1) - (*(pos_j + 1))(1) + ty;
        z[7] = posI2(2) - (*(pos_j + 1))(2) + tz;

        x[8] = posI2(0) - (*(pos_j + 2))(0) + tx;
        y[8] = posI2(1) - (*(pos_j + 2))(1) + ty;
        z[8] = posI2(2) - (*(pos_j + 2))(2) + tz;

        for (int ii = 0; ii < 9; ++ii) {
          r2[ii] = x[ii] * x[ii] + y[ii] * y[ii] + z[ii] * z[ii];
          r2i[ii] = 1.0 / r2[ii];
          ri[ii] = sqrt(r2i[ii]);
        }

        dist6i = r2i[0] * r2i[0] * r2i[0];

        innerloop.spc_innerloop(e_lj, e_crf, dist6i, f, r2, r2i, ri);

        for (int ii = 0; ii < 9; ++ii) {
          fx[ii] = f[ii] * x[ii];
          fy[ii] = f[ii] * y[ii];
          fz[ii] = f[ii] * z[ii];
        }

        (*force_j)(0) -= fx[0] + fx[3] + fx[4];
        (*force_j)(1) -= fy[0] + fy[3] + fy[4];
        (*force_j)(2) -= fz[0] + fz[3] + fz[4];
        (*(force_j + 1))(0) -= fx[1] + fx[5] + fx[7];
        (*(force_j + 1))(1) -= fy[1] + fy[5] + fy[7];
        (*(force_j + 1))(2) -= fz[1] + fz[5] + fz[7];
        (*(force_j + 2))(0) -= fx[2] + fx[6] + fx[8];
        (*(force_j + 2))(1) -= fy[2] + fy[6] + fy[8];
        (*(force_j + 2))(2) -= fz[2] + fz[6] + fz[8];

        groupForce0(0) += fx[0] + fx[1] + fx[2];
        groupForce0(1) += fy[0] + fy[1] + fy[2];
        groupForce0(2) += fz[0] + fz[1] + fz[2];
        groupForce1(0) += fx[3] + fx[5] + fx[6];
        groupForce1(1) += fy[3] + fy[5] + fy[6];
        groupForce1(2) += fz[3] + fz[5] + fz[6];
        groupForce2(0) += fx[4] + fx[7] + fx[8];
        groupForce2(1) += fy[4] + fy[7] + fy[8];
        groupForce2(2) += fz[4] + fz[7] + fz[8];
      }

      storage.force(i) += groupForce0;
      storage.force(i + 1) += groupForce1;
      storage.force(i + 2) += groupForce2;
      
      if (k != 0) {
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
            storage.virial_tensor(b, a) += shift(b) * (groupForce0(a) + groupForce1(a) + groupForce2(a));
          }
        }
      }

      storage.energies.lj_energy[egroup][egroup] += e_lj;
      storage.energies.crf_energy[egroup][egroup] += e_crf;
    }
    
    for (unsigned int ii = end; ii < size_i; ++ii) {
      const math::Vec pos = conf.current().pos(ii);
      const math::Vec force = storage.force(ii);
      //std::cout << "XYZ\t" << ii << std::setprecision(9)
      //          << std::setw(20) << force(0) << std::setw(20) << force(0) << std::setw(20) << force(0)
      //          << std::endl;
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
          storage.virial_tensor(b, a) += pos(b) * force(a);
        }
      }
    }
    
    
  } else if (sim.param().force.special_loop == simulation::special_loop_spc_table) { // special solvent loop
    // solvent - solvent with tabulated spc innerloop...

    // only one energy group
    const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
    if (longrange) {
      double e_lj = 0.0, e_crf = 0.0;
      // here we call the longrange function that uses a smaller table
      for (; i < size_i; i += 3) { // use every third pairlist (OW's)

        const math::Vec posI = conf.current().pos(i);
        const math::Vec posI1 = conf.current().pos(i + 1);
        const math::Vec posI2 = conf.current().pos(i + 2);
        math::Vec groupForce0(0.0);
        math::Vec groupForce1(0.0);
        math::Vec groupForce2(0.0);
        int k = 0;
        math::Vec shift = periodicity.shift(k + 13).pos;
        double tx = shift(0), ty = shift(1), tz = shift(2);

        double r2[9], x[9], y[9], z[9], f[9], fx[9], fy[9], fz[9];
        math::Vec r;

        for (j_it = pairlist_solvent[i].begin(),
            j_to = pairlist_solvent[i].end();
            j_it != j_to;
            j_it += 3) { // use every third atom (OW) in pairlist i

          DEBUG(10, "\tsolvent-solvent (tabulated) longrange spc_nonbonded_interaction: i " << i << " j " << *j_it);

          const int kk = periodicity.nearest_image(posI, conf.current().pos(*j_it), r);

          if (kk != k) {
            storage.force(i) += groupForce0;
            storage.force(i + 1) += groupForce1;
            storage.force(i + 2) += groupForce2;
            if (k != 0) {
              for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                  //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
                  storage.virial_tensor(b, a) += shift(b) * (groupForce0(a) + groupForce1(a) + groupForce2(a));
                }
              }
            }
            groupForce0 = 0.;
            groupForce1 = 0.;
            groupForce2 = 0.;
            k = kk;
            shift = periodicity.shift(k + 13).pos;
            tx = shift(0);
            ty = shift(1);
            tz = shift(2);
          }

          math::Vec const * const pos_j = &conf.current().pos(*j_it);
          math::Vec * const force_j = &storage.force(*j_it);

          x[0] = r(0);
          y[0] = r(1);
          z[0] = r(2);

          x[1] = posI(0) - (*(pos_j + 1))(0) + tx;
          y[1] = posI(1) - (*(pos_j + 1))(1) + ty;
          z[1] = posI(2) - (*(pos_j + 1))(2) + tz;

          x[2] = posI(0) - (*(pos_j + 2))(0) + tx;
          y[2] = posI(1) - (*(pos_j + 2))(1) + ty;
          z[2] = posI(2) - (*(pos_j + 2))(2) + tz;

          x[3] = posI1(0) - (*(pos_j))(0) + tx;
          y[3] = posI1(1) - (*(pos_j))(1) + ty;
          z[3] = posI1(2) - (*(pos_j))(2) + tz;

          x[4] = posI2(0) - (*(pos_j))(0) + tx;
          y[4] = posI2(1) - (*(pos_j))(1) + ty;
          z[4] = posI2(2) - (*(pos_j))(2) + tz;

          x[5] = posI1(0) - (*(pos_j + 1))(0) + tx;
          y[5] = posI1(1) - (*(pos_j + 1))(1) + ty;
          z[5] = posI1(2) - (*(pos_j + 1))(2) + tz;

          x[6] = posI1(0) - (*(pos_j + 2))(0) + tx;
          y[6] = posI1(1) - (*(pos_j + 2))(1) + ty;
          z[6] = posI1(2) - (*(pos_j + 2))(2) + tz;

          x[7] = posI2(0) - (*(pos_j + 1))(0) + tx;
          y[7] = posI2(1) - (*(pos_j + 1))(1) + ty;
          z[7] = posI2(2) - (*(pos_j + 1))(2) + tz;

          x[8] = posI2(0) - (*(pos_j + 2))(0) + tx;
          y[8] = posI2(1) - (*(pos_j + 2))(1) + ty;
          z[8] = posI2(2) - (*(pos_j + 2))(2) + tz;

          for (int ii = 0; ii < 9; ++ii) {
            r2[ii] = x[ii] * x[ii] + y[ii] * y[ii] + z[ii] * z[ii];
          }

          innerloop.longrange_spc_table_innerloop(e_lj, e_crf, f, r2);

          for (int ii = 0; ii < 9; ++ii) {
            fx[ii] = f[ii] * x[ii];
            fy[ii] = f[ii] * y[ii];
            fz[ii] = f[ii] * z[ii];
          }

          (*force_j)(0) -= fx[0] + fx[3] + fx[4];
          (*force_j)(1) -= fy[0] + fy[3] + fy[4];
          (*force_j)(2) -= fz[0] + fz[3] + fz[4];
          (*(force_j + 1))(0) -= fx[1] + fx[5] + fx[7];
          (*(force_j + 1))(1) -= fy[1] + fy[5] + fy[7];
          (*(force_j + 1))(2) -= fz[1] + fz[5] + fz[7];
          (*(force_j + 2))(0) -= fx[2] + fx[6] + fx[8];
          (*(force_j + 2))(1) -= fy[2] + fy[6] + fy[8];
          (*(force_j + 2))(2) -= fz[2] + fz[6] + fz[8];

          groupForce0(0) += fx[0] + fx[1] + fx[2];
          groupForce0(1) += fy[0] + fy[1] + fy[2];
          groupForce0(2) += fz[0] + fz[1] + fz[2];
          groupForce1(0) += fx[3] + fx[5] + fx[6];
          groupForce1(1) += fy[3] + fy[5] + fy[6];
          groupForce1(2) += fz[3] + fz[5] + fz[6];
          groupForce2(0) += fx[4] + fx[7] + fx[8];
          groupForce2(1) += fy[4] + fy[7] + fy[8];
          groupForce2(2) += fz[4] + fz[7] + fz[8];
        }

        storage.force(i) += groupForce0;
        storage.force(i + 1) += groupForce1;
        storage.force(i + 2) += groupForce2;
        if (k != 0) {
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
              storage.virial_tensor(b, a) += shift(b) * (groupForce0(a) + groupForce1(a) + groupForce2(a));
            }
          }
        }
      }

    for (unsigned int ii = topo.num_solute_atoms(); ii < size_i; ++ii) {
      const math::Vec pos = conf.current().pos(ii);
      const math::Vec force = storage.force(ii);
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          //innerloop.spc_innerloop(topo, conf, i, *j_it, storage, periodicity);
          storage.virial_tensor(b, a) += pos(b) * force(a);
        }
      }
    }

      storage.energies.lj_energy[egroup][egroup] += e_lj;
      storage.energies.crf_energy[egroup][egroup] += e_crf;
    } else { // shortrange
      double e_lj = 0.0, e_crf = 0.0;
      for (; i < size_i; i += 3) { // use every third pairlist (OW's)

        const math::Vec posI = conf.current().pos(i);
        const math::Vec posI1 = conf.current().pos(i + 1);
        const math::Vec posI2 = conf.current().pos(i + 2);
        math::Vec groupForce0(0.0);
        math::Vec groupForce1(0.0);
        math::Vec groupForce2(0.0);
        int k = 0;
        math::Vec shift = periodicity.shift(k + 13).pos;
        double tx = shift(0), ty = shift(1), tz = shift(2);

        double r2[9], x[9], y[9], z[9], f[9], fx[9], fy[9], fz[9];
        math::Vec r;

        for (j_it = pairlist_solvent[i].begin(),
                j_to = pairlist_solvent[i].end();
                j_it != j_to;
                j_it += 3) { // use every third atom (OW) in pairlist i

          DEBUG(10, "\tsolvent-solvent (tabulated) shortrange spc_nonbonded_interaction: i " << i << " j " << *j_it);

          const int kk = periodicity.nearest_image(posI, conf.current().pos(*j_it), r);

          if (kk != k) {
            storage.force(i) += groupForce0;
            storage.force(i + 1) += groupForce1;
            storage.force(i + 2) += groupForce2;
            if (k != 0) {
              for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                  //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
                  storage.virial_tensor(b, a) += shift(b) * (groupForce0(a) + groupForce1(a) + groupForce2(a));
                }
              }
            }
            groupForce0 = 0.;
            groupForce1 = 0.;
            groupForce2 = 0.;
            k = kk;
            shift = periodicity.shift(k + 13).pos;
            tx = shift(0);
            ty = shift(1);
            tz = shift(2);
          }

          math::Vec const * const pos_j = &conf.current().pos(*j_it);
          math::Vec * const force_j = &storage.force(*j_it);

          x[0] = r(0);
          y[0] = r(1);
          z[0] = r(2);

          x[1] = posI(0) - (*(pos_j + 1))(0) + tx;
          y[1] = posI(1) - (*(pos_j + 1))(1) + ty;
          z[1] = posI(2) - (*(pos_j + 1))(2) + tz;

          x[2] = posI(0) - (*(pos_j + 2))(0) + tx;
          y[2] = posI(1) - (*(pos_j + 2))(1) + ty;
          z[2] = posI(2) - (*(pos_j + 2))(2) + tz;

          x[3] = posI1(0) - (*(pos_j))(0) + tx;
          y[3] = posI1(1) - (*(pos_j))(1) + ty;
          z[3] = posI1(2) - (*(pos_j))(2) + tz;

          x[4] = posI2(0) - (*(pos_j))(0) + tx;
          y[4] = posI2(1) - (*(pos_j))(1) + ty;
          z[4] = posI2(2) - (*(pos_j))(2) + tz;

          x[5] = posI1(0) - (*(pos_j + 1))(0) + tx;
          y[5] = posI1(1) - (*(pos_j + 1))(1) + ty;
          z[5] = posI1(2) - (*(pos_j + 1))(2) + tz;

          x[6] = posI1(0) - (*(pos_j + 2))(0) + tx;
          y[6] = posI1(1) - (*(pos_j + 2))(1) + ty;
          z[6] = posI1(2) - (*(pos_j + 2))(2) + tz;

          x[7] = posI2(0) - (*(pos_j + 1))(0) + tx;
          y[7] = posI2(1) - (*(pos_j + 1))(1) + ty;
          z[7] = posI2(2) - (*(pos_j + 1))(2) + tz;

          x[8] = posI2(0) - (*(pos_j + 2))(0) + tx;
          y[8] = posI2(1) - (*(pos_j + 2))(1) + ty;
          z[8] = posI2(2) - (*(pos_j + 2))(2) + tz;

          for (int ii = 0; ii < 9; ++ii) {
            r2[ii] = x[ii] * x[ii] + y[ii] * y[ii] + z[ii] * z[ii];
          }

          innerloop.shortrange_spc_table_innerloop(e_lj, e_crf, f, r2);

          for (int ii = 0; ii < 9; ++ii) {
            fx[ii] = f[ii] * x[ii];
            fy[ii] = f[ii] * y[ii];
            fz[ii] = f[ii] * z[ii];
          }

          (*force_j)(0) -= fx[0] + fx[3] + fx[4];
          (*force_j)(1) -= fy[0] + fy[3] + fy[4];
          (*force_j)(2) -= fz[0] + fz[3] + fz[4];
          (*(force_j + 1))(0) -= fx[1] + fx[5] + fx[7];
          (*(force_j + 1))(1) -= fy[1] + fy[5] + fy[7];
          (*(force_j + 1))(2) -= fz[1] + fz[5] + fz[7];
          (*(force_j + 2))(0) -= fx[2] + fx[6] + fx[8];
          (*(force_j + 2))(1) -= fy[2] + fy[6] + fy[8];
          (*(force_j + 2))(2) -= fz[2] + fz[6] + fz[8];

          groupForce0(0) += fx[0] + fx[1] + fx[2];
          groupForce0(1) += fy[0] + fy[1] + fy[2];
          groupForce0(2) += fz[0] + fz[1] + fz[2];
          groupForce1(0) += fx[3] + fx[5] + fx[6];
          groupForce1(1) += fy[3] + fy[5] + fy[6];
          groupForce1(2) += fz[3] + fz[5] + fz[6];
          groupForce2(0) += fx[4] + fx[7] + fx[8];
          groupForce2(1) += fy[4] + fy[7] + fy[8];
          groupForce2(2) += fz[4] + fz[7] + fz[8];
        }

        storage.force(i) += groupForce0;
        storage.force(i + 1) += groupForce1;
        storage.force(i + 2) += groupForce2;
        if (k != 0) {
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
              storage.virial_tensor(b, a) += shift(b) * (groupForce0(a) + groupForce1(a) + groupForce2(a));
            }
          }
        }
      }

    for (unsigned int ii = topo.num_solute_atoms(); ii < size_i; ++ii) {
      const math::Vec pos = conf.current().pos(ii);
      const math::Vec force = storage.force(ii);
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          //innerloop.spc_innerloop(topo, conf, i, *j_it, storage, periodicity);
          storage.virial_tensor(b, a) += pos(b) * force(a);
        }
      }
    }

      storage.energies.lj_energy[egroup][egroup]  += e_lj;
      storage.energies.crf_energy[egroup][egroup] += e_crf;
    }
  } else if (sim.param().force.special_loop == simulation::special_loop_generic) {
    const unsigned int num_solvent_atoms = topo.solvent(0).num_atoms();
    // prepare parameters
    const unsigned int num_param = num_solvent_atoms * num_solvent_atoms;
    typename interaction::Nonbonded_Innerloop<interaction::Interaction_Spec<math::rectangular, simulation::lj_crf_func> >::solvent_pair_parameter pair_parameter[num_param];

    unsigned int param = 0;
    for (unsigned int atom_i = 0; atom_i < num_solvent_atoms; ++atom_i) {
      for (unsigned int atom_j = 0; atom_j < num_solvent_atoms; ++atom_j, ++param) {
        assert(param < num_param);
        DEBUG(10, "\tsolvent pair parameter: " << param);

        const lj_parameter_struct & lj =
                innerloop.param()->lj_parameter(topo.solvent(0).atom(atom_i).iac, topo.solvent(0).atom(atom_j).iac);
        pair_parameter[param].c12 = lj.c12;
        pair_parameter[param].c6 = lj.c6;

        pair_parameter[param].q = math::four_pi_eps_i *
                topo.solvent(0).atom(atom_i).charge *
                topo.solvent(0).atom(atom_j).charge;

        DEBUG(10, "\t\tc12: " << pair_parameter[param].c12 << " c6: " <<
                pair_parameter[param].c6 << " q: " << pair_parameter[param].q);
      }
    }

    math::Vec r;
    
    const double crf_2cut3i = innerloop.crf_2cut3i(0);
    const double crf_cut3i  = crf_2cut3i * 2.;
    const double crf_cut    = innerloop.crf_cut(0);
  
    // only one energy group
    const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
    double e_lj = 0.0, e_crf = 0.0;
    // use num_solvent_atoms-th atom (first of solvent molecule i)
    for (; i < size_i; i += num_solvent_atoms) {
      math::Vec const * const pos_i = &conf.current().pos(i);
      math::Vec * const force_i = &storage.force(i);
      int k = 0;
      std::vector<math::Vec> groupForce(num_solvent_atoms, math::Vec(0.0));
      math::Vec totalForce;
      math::Vec shift = periodicity.shift(k+13).pos;
      double tx = shift(0), ty = shift(1), tz = shift(2);
      
      for (j_it = pairlist_solvent[i].begin(),
              j_to = pairlist_solvent[i].end();
              j_it != j_to;
              j_it += num_solvent_atoms) { // use num_solvent_atoms-th atom (first of solvent molecule j)

        DEBUG(10, "\tsolvent_nonbonded_interaction (special_loop_generic): i " << i << " j " << *j_it);

        const int kk = periodicity.nearest_image(*pos_i, conf.current().pos(*j_it), r);
        
        if (kk != k) {
	  totalForce = 0.;
          for(unsigned int atom_i = 0; atom_i < num_solvent_atoms; ++atom_i) {
            (*(force_i+atom_i)) += groupForce[atom_i];
            if (k != 0) {
	      totalForce += groupForce[atom_i];
            }            
            groupForce[atom_i] = 0.0;
          }
	  if (k != 0) {
            const math::Vec shift = periodicity.shift(k+13).pos;
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) {
                //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
                storage.virial_tensor(b, a) += shift(b) * totalForce(a);
              }
            }
	  }
          k = kk;
          shift = periodicity.shift(k + 13).pos;
          tx = shift(0);
          ty = shift(1);
          tz = shift(2);
        }
        
        math::Vec const * const pos_j = &conf.current().pos(*j_it);
        math::Vec * const force_j = &storage.force(*j_it);
        
        for(unsigned int param = 0, atom_i = 0; atom_i < num_solvent_atoms; ++atom_i) {
          const double xi = (*(pos_i+atom_i))(0) + tx;
          const double yi = (*(pos_i+atom_i))(1) + ty;
          const double zi = (*(pos_i+atom_i))(2) + tz;
          for(unsigned int atom_j = 0; atom_j < num_solvent_atoms; ++atom_j, ++param) {
            DEBUG(15, "\tatoms: i: " << atom_i << " j: " << atom_j);
            const double x = xi - (*(pos_j+atom_j))(0);
            const double y = yi - (*(pos_j+atom_j))(1);
            const double z = zi - (*(pos_j+atom_j))(2);

            const double r2 = x * x + y * y + z * z;
            DEBUG(15, "\tr2: " << r2);
            assert(r2 != 0.0);
            const double r2i = 1.0 / r2;
            const double ri = sqrt(r2i);
            const double dist6i = r2i * r2i * r2i;
            const double dist6i_c12 = pair_parameter[param].c12 * dist6i;

            e_lj += (dist6i_c12 - pair_parameter[param].c6) * dist6i;
            e_crf += pair_parameter[param].q * (ri - crf_2cut3i * r2 - crf_cut);

            const double f = (dist6i_c12 + dist6i_c12 - pair_parameter[param].c6) * 6.0 
                    * dist6i * r2i + pair_parameter[param].q * (ri * r2i + crf_cut3i);

            const double fx = f * x;
            const double fy = f * y;
            const double fz = f * z;

            //(*(force_i+atom_i))(0) += fx;
            groupForce[atom_i](0) += fx;
            (*(force_j+atom_j))(0) -= fx;
            //(*(force_i+atom_i))(1) += fy;
            groupForce[atom_i](1) += fy;
            (*(force_j+atom_j))(1) -= fy;
            //(*(force_i+atom_i))(2) += fz;
            groupForce[atom_i](2) += fz;
            (*(force_j+atom_j))(2) -= fz;
          }
        // innerloop.solvent_innerloop(topo, pair_parameter, conf, num_solvent_atoms, i, *j_it, storage, periodicity);
        }
      }
        
      totalForce = 0.;
      for(unsigned int atom_i = 0; atom_i < num_solvent_atoms; ++atom_i) {
        (*(force_i+atom_i)) += groupForce[atom_i];
        if (k != 0) {
	  totalForce += groupForce[atom_i];
        }
      }
        
      if (k != 0) {
        const math::Vec shift = periodicity.shift(k+13).pos;
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
            storage.virial_tensor(b, a) += shift(b) * totalForce(a);
          }
        }
      }   
    }
    
    storage.energies.lj_energy[egroup][egroup] += e_lj;
    storage.energies.crf_energy[egroup][egroup] += e_crf;
    
    
    for (unsigned int ii = end; ii < size_i; ++ii) {
      const math::Vec pos = conf.current().pos(ii);
      const math::Vec force = storage.force(ii);
      //std::cout << "XYZ\t" << ii << std::setprecision(9)
      //          << std::setw(20) << force(0) << std::setw(20) << force(0) << std::setw(20) << force(0)
      //          << std::endl;
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
          storage.virial_tensor(b, a) += pos(b) * force(a);
        }
      }
    }
    
    
  } else { // normal solvent loop
    for (; i < size_i; ++i) {
      const math::Vec posI = conf.current().pos(i);
      const unsigned int eg_i = topo.atom_energy_group(i);
      math::Vec groupForce(0.0);
      int k = 0;
      for (j_it = pairlist_solvent[i].begin(),
              j_to = pairlist_solvent[i].end();
              j_it != j_to;
              ++j_it) {

        DEBUG(10, "\tsolvent_nonbonded_interaction (normal solvent loop): i " << i << " j " << *j_it);
            
        math::Vec r;
        const int kk = periodicity.nearest_image(posI, conf.current().pos(*j_it), r);

        if (kk != k) {
          storage.force(i) += groupForce;
          if (k != 0) {
            const math::Vec shift = periodicity.shift(k+13).pos;
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) {
                //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
                storage.virial_tensor(b, a) += shift(b) * groupForce(a);
              }
            }
          }
          groupForce = 0.;
          k = kk;

        }

        const double dist2 = abs2(r);
        math::Vec force;
        double f = 0.0;
        double e_lj = 0.0, e_crf = 0.0;

        innerloop.lj_crf_innerloop_2(topo, i, *j_it, dist2, f, e_lj, e_crf);

        const unsigned int eg_j = topo.atom_energy_group(*j_it);
        DEBUG(11, "\tenergy group i " << eg_i << " j " << eg_j);
        storage.energies.lj_energy[eg_i][eg_j] += e_lj;
        storage.energies.crf_energy[eg_i][eg_j] += e_crf;

        force = f * r;
        storage.force(*j_it) -= force;
        groupForce += force;

        //innerloop.lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
      }
    
      storage.force(i) += groupForce;
      if (k != 0) {
        const math::Vec shift = periodicity.shift(k+13).pos;
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) {
            //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
            storage.virial_tensor(b, a) += shift(b) * groupForce(a);
          }
        }
      }
    }

    for (unsigned int ii = end; ii < size_i; ++ii) {
      const math::Vec pos = conf.current().pos(ii);
      const math::Vec force = storage.force(ii);
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          //storage.virial_tensor(b, a) += (posI(b) + shift(b)) * groupForce(a);
          storage.virial_tensor(b, a) += pos(b) * force(a);
        }
      }
    }
  }
  if (master)
    timer.stop(timer_name);
}

void interaction::Nonbonded_Outerloop
::lj_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Pairlist const & pairlist_solute,
        Pairlist const & pairlist_solvent,
        Storage & storage,
        bool longrange, util::Algorithm_Timer & timer, bool master) {
  SPLIT_MY_INNERLOOP(_lj_outerloop, simulation::lj_func, topo, conf, sim,
          pairlist_solute, pairlist_solvent, storage, longrange, timer, master);
}

template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_lj_outerloop(topology::Topology & topo,
                configuration::Configuration & conf,
                simulation::Simulation & sim,
                Pairlist const & pairlist_solute,
                Pairlist const & pairlist_solvent,
                Storage & storage,
                bool longrange, util::Algorithm_Timer & timer, bool master) {
  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);

  innerloop.init(sim,simulation::lj_func);

  const unsigned end_solute = topo.num_solute_atoms();
  for (unsigned i = 0; i < end_solute; ++i) {
    for (std::vector<unsigned int>::const_iterator
          j_it = pairlist_solute[i].begin()
        , j_to = pairlist_solute[i].end()
        ; j_it != j_to; ++j_it) {
      DEBUG(10, "\tLJ nonbonded_interaction: i " << i << " j " << *j_it);

      innerloop.lj_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }
  const unsigned end_solvent = topo.num_atoms();
  for (unsigned i = end_solute; i < end_solvent; ++i) {
    for (std::vector<unsigned int>::const_iterator
          j_it = pairlist_solvent[i].begin()
        , j_to = pairlist_solvent[i].end()
        ; j_it != j_to; ++j_it) {
      DEBUG(10, "\tLJ nonbonded_interaction: i " << i << " j " << *j_it);

      innerloop.lj_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }
}

// calculate sasa and volume term

void interaction::Nonbonded_Outerloop
::sasa_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage) {
  SPLIT_INNERLOOP(_sasa_outerloop, topo, conf, sim,
          storage);
}

/**
 * helper function to calculate sasa forces and energies
 */

template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_sasa_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage) {

  DEBUG(7, "\tCalculating SASA/VOL interaction term");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);

  // these are the same for every atom, but I only want them if it's volume...
  const double amin = sim.param().sasa.min_cut;
  const double amax = sim.param().sasa.max_cut;
  const double adiff = sim.param().sasa.cut_diff;

  // first put surface and area into conf.current
  const unsigned int num_sasa_atoms = topo.sasa_parameter().size();
  DEBUG(15, "\tNumber of non-H (\"sasa\") atoms: " << num_sasa_atoms);

  for (unsigned int i = 0; i < num_sasa_atoms; ++i) {

    // note that i counts non-H ("sasa") atoms
    DEBUG(10, "\tInitialising surface, area and volume for sasa atom " << i);
    // get sasa parameters for atom i
    const topology::sasa_parameter_struct & sasa_param_i = topo.sasa_parameter(i);
    // initialise actual sasa array (size = num_sasa_atoms)
    conf.current().sasa_area[i] = sasa_param_i.surface;

  } // end initialize surface areas

  // now compute actual sasa (reduction due to overlap)
  for (unsigned int i = 0; i < num_sasa_atoms; ++i) {

    DEBUG(10, "\tCalculating true SASA for sasa atom " << i);
    innerloop.sasa_calc_innerloop(topo, conf, i, sim, periodicity);

  } // end compute sasa

  // store final and total sasas and their energy contribution,
  // compute volume term and compute forces
  for (unsigned int i = 0; i < num_sasa_atoms; ++i) {

    // check for negative SASA
    DEBUG(10, "\tChecking for negative true SASA for sasa atom " << i);
    if (conf.current().sasa_area[i] < 0) {
      io::messages.add("Nonbonded_Outerloop",
              "negative SASA", io::message::critical);
    } else {
      DEBUG(10, "\tStoring true SASA and energies for sasa atom " << i);

      // get sasa parameters for atom i
      const topology::sasa_parameter_struct & sasa_param_i = topo.sasa_parameter(i);
      // add sasa to total
      conf.current().sasa_tot += conf.current().sasa_area[i];

      DEBUG(15, "\tSASA of atom " << sasa_param_i.atom << " is " << conf.current().sasa_area[i]
              << "\tand current total SASA of molecule is " << conf.current().sasa_tot);

      double e_sasa = conf.current().sasa_area[i] * sasa_param_i.sigma;
      conf.current().energies.sasa_energy[topo.atom_energy_group(sasa_param_i.atom)] += e_sasa;

      DEBUG(15, "\tSASA energy of atom " << sasa_param_i.atom << " is " << e_sasa
              << "\tand current total SASA energy is " <<
              conf.current().energies.sasa_energy[topo.atom_energy_group(sasa_param_i.atom)]);

      // if using volume too, compute volume term for atom i
      // has to go here because the switching function needs the true area
      if (sim.param().sasa.switch_volume) {
        DEBUG(10, "\tComputing volume term for sasa atom " << i);
        // first compute switching function and its derivative (for later)
        innerloop.sasa_switching_fct(conf, i, amin, amax, adiff);
        // compute volume contribution to energy
        innerloop.sasa_volume_innerloop(topo, conf, i, sim);
      }
    } // end else
  } // end atoms i

  // finally calculate the forces (for sasa and, if used, vol)
  // has to be in a separate loop because we need the switching function derivative
  // for all atoms to have been computed
  for (unsigned int i = 0; i < num_sasa_atoms; ++i) {
    DEBUG(10, "\tCalculating SASA/VOL forces for sasa atom " << i);
    innerloop.sasa_force_innerloop(topo, conf, i, sim, periodicity);

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
        Storage & storage,
        int rank, int size) {
  SPLIT_INNERLOOP(_one_four_outerloop, topo, conf, sim, storage, rank, size);
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
        Storage & storage,
        int rank, int size) {
  DEBUG(7, "\tcalculate 1,4-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  topology::excl_cont_t::value_type::const_iterator it, to;
  unsigned int const num_solute_atoms = topo.num_solute_atoms();

  for (unsigned int i = rank; i < num_solute_atoms; i += size) {
    it = topo.one_four_pair(i).begin();
    to = topo.one_four_pair(i).end();

    for (; it != to; ++it) {

      innerloop.one_four_interaction_innerloop(topo, conf, i, *it, storage, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}

/**
 * helper function to calculate the forces and energies from the
 * Lennard-Jones exception interaction
 */
void interaction::Nonbonded_Outerloop
::lj_exception_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage,
                     int rank, int size)
{
  SPLIT_INNERLOOP(_lj_exception_outerloop, topo, conf, sim, storage, rank, size);
}

/**
 * helper function to calculate the forces and energies from the
 * Lennard-Jones exception interaction
 */
template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_lj_exception_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage,
                     int rank, int size)
{
  DEBUG(7, "\tcalculate Lennard-Jones-exception-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  unsigned int const num_lj_exceptions = topo.lj_exceptions().size();

  for (unsigned int i = rank; i < num_lj_exceptions; i += size) {
    const topology::lj_exception_struct & ljex = topo.lj_exceptions()[i];

    innerloop.lj_exception_innerloop(topo, conf, ljex, storage, periodicity);
  } // loop over LJ exceptions
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
void interaction::Nonbonded_Outerloop
::cg_exclusions_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage) {
  SPLIT_INNERLOOP(_cg_exclusions_outerloop, topo, conf, sim, storage);
}

template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_cg_exclusions_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage) {
  DEBUG(7, "\tcalculate 1,4-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  std::set<int>::const_iterator cg2_it, cg2_to;

  for (unsigned int cg1 = 0; cg1 < topo.num_solute_chargegroups(); ++cg1) {

    cg2_it = topo.chargegroup_exclusion(cg1).begin();
    cg2_to = topo.chargegroup_exclusion(cg1).end();

    for (; cg2_it != cg2_to; ++cg2_it) {

      // only once...
      if (cg1 > (unsigned) * cg2_it) continue;

      for (int a1 = topo.chargegroup(cg1); a1 < topo.chargegroup(cg1 + 1); ++a1) {
        for (int a2 = topo.chargegroup(*cg2_it); a2 < topo.chargegroup(*cg2_it + 1); ++a2) {

          // std::cout << "cg1=" << cg1 << " cg2=" << *cg2_it
          // << " a1=" << a1 << " a2=" << a2 << std::endl;

          if (a1 >= a2) continue;

          if (std::find(topo.exclusion(a1).begin(), topo.exclusion(a1).end(), a2) 
              != topo.exclusion(a1).end()) continue;

          if (std::find(topo.one_four_pair(a1).begin(), topo.one_four_pair(a1).end(), a2) 
              != topo.one_four_pair(a1).end()) {
            // std::cout << "\t1,4" << std::endl;
            innerloop.one_four_interaction_innerloop(topo, conf, a1, a2, storage, periodicity);
          } else {
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
        Storage & storage,
        int rank, int size) {
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

  SPLIT_INNERLOOP(_RF_excluded_outerloop, topo, conf, sim, storage, rank, size);
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
        Storage & storage,
        int rank, int size) {

  DEBUG(7, "\tcalculate RF excluded interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  const unsigned int num_solute_atoms = topo.num_solute_atoms();

  // Skip QM atoms in electrostatic or polarisable embedding
  const bool qmmm = (sim.param().qmmm.qmmm > simulation::qmmm_mechanical);
  for (unsigned int i = rank; i < num_solute_atoms; i += size) {
    if (qmmm && topo.is_qm(i)) continue;
    innerloop.RF_excluded_interaction_innerloop(topo, conf, i, storage, periodicity);
  } // loop over solute atoms

  // Solvent
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
          cg_to = topo.chargegroup_end();
  cg_it += topo.num_solute_chargegroups() + rank;

  for (; cg_it < cg_to; cg_it += size) {
    innerloop.RF_solvent_interaction_innerloop(topo, conf, cg_it, storage, periodicity);
  } // loop over solvent charge groups
}

void interaction::Nonbonded_Outerloop
::self_energy_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage) {
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
        Storage & storage) {
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
    if (topo.is_polarisable(i)) {

      DEBUG(10, "\tself energy: i " << i);
      innerloop.self_energy_innerloop(
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
        Storage & storage_lr, int rank) {
  SPLIT_INNERLOOP(_electric_field_outerloop, topo, conf, sim,
          pairlist, storage, storage_lr, rank);
}

/**
 * helper function to calculate polarisation, 
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
        int rank) {
  DEBUG(7, "\tcalculate polarisation (electric field outerloop)");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  unsigned int i = 0;
  unsigned int size_i = unsigned(pairlist.size());
  unsigned int size_lr = size_i;
  DEBUG(11, "el_field outerloop pairlist size " << size_i);

  unsigned int end = size_i;
  unsigned int end_lr = size_lr;
  
  const bool do_qmmm = 
      (rank == 0 && sim.param().qmmm.qmmm == simulation::qmmm_polarisable);

  QMMM_Interaction * qmmm = nullptr;
  if (do_qmmm) {
    qmmm = QMMM_Interaction::pointer();
    if (qmmm == nullptr) {
      io::messages.add("Unable to get QMMM interaction in electric field calculation",
                        "Nonbonded_Outerloop", io::message::critical);
    }
  }

  math::VArray e_el_new(topo.num_atoms());
#ifdef OMP
#pragma omp single
  {
    Nonbonded_Outerloop::electric_field.resize(e_el_new.size());
    Nonbonded_Outerloop::electric_field = 0.0;
  }
#pragma omp barrier
#endif

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
  if (!(sim.steps() % sim.param().pairlist.skip_step)) {

    // loop over all molecules in longrange pairlist
    for (i = 0; i < end_lr; ++i) {
      // solute longrange
      for (std::vector<unsigned int>::const_iterator
        j_it = pairlist.solute_long[i].begin(),
              j_to = pairlist.solute_long[i].end();
              j_it != j_to; ++j_it){

        if (!topo.is_polarisable(i) && !topo.is_polarisable(*j_it)) continue;
        math::Vec e_eli_lr, e_elj_lr;

        innerloop.electric_field_innerloop(topo, conf, i, *j_it,
                e_eli_lr, e_elj_lr, periodicity);

        storage_lr.electric_field[i] += e_eli_lr;
        storage_lr.electric_field[*j_it] += e_elj_lr;
      }
      // solvent longrange
      for (std::vector<unsigned int>::const_iterator
        j_it = pairlist.solvent_long[i].begin(),
              j_to = pairlist.solvent_long[i].end();
              j_it != j_to; ++j_it){
        if (!topo.is_polarisable(i) && !topo.is_polarisable(*j_it)) continue;
        math::Vec e_eli_lr, e_elj_lr;

        innerloop.electric_field_innerloop(topo, conf, i, *j_it,
                e_eli_lr, e_elj_lr, periodicity);

        storage_lr.electric_field[i] += e_eli_lr;
        storage_lr.electric_field[*j_it] += e_elj_lr;
      }
    }

// Reduce longrange
#ifdef OMP
    // Reduce longrange electric field from OMP threads
    for (unsigned i = 0; i < topo.num_atoms(); ++i) {
      for (unsigned j = 0; j < 3; ++j) {
#pragma omp atomic
        Nonbonded_Outerloop::electric_field[i][j] += storage_lr.electric_field[i][j];
      }
    }
// Wait for all writes
#pragma omp barrier
    if (rank == 0) {
      storage_lr.electric_field.swap(Nonbonded_Outerloop::electric_field);
      Nonbonded_Outerloop::electric_field = 0.0;
    }
#pragma omp barrier
#endif

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
    for (i = 0; i < end; ++i) {
      // solute short
      for (std::vector<unsigned int>::const_iterator
        j_it = pairlist.solute_short[i].begin(),
              j_to = pairlist.solute_short[i].end();
              j_it != j_to; ++j_it) {
        if (!topo.is_polarisable(i) && !topo.is_polarisable(*j_it)) continue;
        math::Vec e_eli, e_elj;

        innerloop.electric_field_innerloop(topo, conf,
                i, *j_it, e_eli, e_elj, periodicity);

        e_el_new(i) += e_eli;
        e_el_new(*j_it) += e_elj;
      }
      // solvent short
      for (std::vector<unsigned int>::const_iterator
        j_it = pairlist.solvent_short[i].begin(),
              j_to = pairlist.solvent_short[i].end();
              j_it != j_to; ++j_it) {
        if (!topo.is_polarisable(i) && !topo.is_polarisable(*j_it)) continue;
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

#ifdef OMP
    // Reduce electric field from OMP threads
    for (unsigned i = 0; i < topo.num_atoms(); ++i) {
      for (unsigned j = 0; j < 3; ++j) {
#pragma omp atomic
        Nonbonded_Outerloop::electric_field[i][j] += e_el_new[i][j];
      }
    }
// Wait for all writes
#pragma omp barrier
    if (rank == 0) {
      e_el_new.swap(Nonbonded_Outerloop::electric_field);
      Nonbonded_Outerloop::electric_field = 0.0;
    }
#endif

    if (rank == 0) {
      if (do_qmmm) {
        // get the contributions from the QM part
        if (turni > 0) {
          // First QM iteration was already performed
          qmmm->scf_step(topo, conf, sim);
        }
        qmmm->get_electric_field(sim, e_el_new);
      }
      
      // If the external electric field is activated
      // then it will also act on the polarisable model
      // as a perturbation to the COS electric field
      math::Vec external_field(sim.param().electric.Ef_x,
                               sim.param().electric.Ef_y,
                               sim.param().electric.Ef_z);
      external_field *= math::four_pi_eps_i * 3 * sim.param().nonbonded.rf_epsilon / (2 * sim.param().nonbonded.rf_epsilon
              + sim.param().nonbonded.epsilon);
      
      for (i=0; i<topo.num_atoms(); ++i) {
        if(topo.is_polarisable(i)){
          e_el_new(i) += storage_lr.electric_field(i) + (external_field);
          
          //delta r
          math::Vec delta_r;

          //////////////////////////////////////////////////
          // implementation of polarisability damping
          /////////////////////////////////////////////////

          if (sim.param().polarise.damp) { // damp the polarisability
            const double e_i = math::abs(e_el_new(i)),
                    e_0 = topo.damping_level(i);
            if (e_i <= e_0)
              delta_r = (topo.polarisability(i) / topo.coscharge(i)) * e_el_new(i);
            else {
              const double p = topo.damping_power(i);
              delta_r = topo.polarisability(i) * e_0  * 
                        (p + 1.0 - pow(e_0/e_i, p)) / 
                        (p * topo.coscharge(i) * e_i) * e_el_new(i);
            }
          } else { // no damping
            delta_r = (topo.polarisability(i) / topo.coscharge(i)) * e_el_new(i);
          }
          // store the new position
          conf.current().posV(i) = delta_r;

          // calculation of convergence criterium
          for (int j = 0; j < 3; ++j) {
            double delta_field = fabs(storage.electric_field(i)(j) - e_el_new(i)(j));
            if (delta_field > maxfield) {
              maxfield = delta_field;
            }
          }
        }

        storage.electric_field(i) = e_el_new(i);
      }
    } // end if rank==0
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
#ifdef OMP
  // share the convergence criterion
    if (rank == 0)
      Nonbonded_Outerloop::minfield = minfield;
#pragma omp barrier
    if (rank != 0)
      minfield = Nonbonded_Outerloop::minfield;
#endif

    DEBUG(11, "\trank: " << rank << " minfield: " << minfield << " iteration round: " << turni);
  }

  if (do_qmmm)
    qmmm->write_qm_data(topo, conf, sim);
  DEBUG(5, "electric field iterations: " << turni);
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
        Storage & storage, int rank, int size) {

  DEBUG(7, "\tcalculate LS real space interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  /*
    variables for a OMP parallelizable loop.
    outer index has to be integer...
   */
  std::vector<unsigned int>::const_iterator j_it, j_to;
  topology::excl_cont_t::value_type::const_iterator ex_it, ex_to;

  unsigned int size_i = unsigned(pairlist_solute.size());
  DEBUG(10, "ls_real outerloop pairlist size " << size_i);

  unsigned int end = topo.num_solute_atoms();

  unsigned int i = 0;
  for (i = 0; i < end; i++) {
    for (j_it = pairlist_solute[i].begin(),
            j_to = pairlist_solute[i].end();
            j_it != j_to;
            ++j_it) {

      DEBUG(10, "\tnonbonded_interaction: i " << i << " j " << *j_it);

      // shortrange, therefore store in simulation.system()
      innerloop.lj_ls_real_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }

  for (; i < size_i; i++) {
    for (j_it = pairlist_solvent[i].begin(),
            j_to = pairlist_solvent[i].end();
            j_it != j_to;
            ++j_it) {

      DEBUG(10, "\tsolvent_nonbonded_interaction: i " << i << " j " << *j_it);

      innerloop.lj_ls_real_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }

  // loop over all exclusions as they have reduced LS interactions in real space
  // parrallelization using stride. Then MPI should work
  DEBUG(9, "U_eta due to excluded solute pairs...");
  const unsigned int size_int = topo.num_solute_atoms();
  for (unsigned int i = rank; i < size_int; i += size) {
    for (ex_it = topo.exclusion(i).begin(),
            ex_to = topo.exclusion(i).end();
            ex_it != ex_to;
            ++ex_it) {

      DEBUG(10, "\texcluded_nonbonded_interaction: i " << i << " j " << *ex_it);
      innerloop.ls_real_excluded_innerloop(topo, conf, i, *ex_it, storage, periodicity);
    }
  }
  DEBUG(9, "U_eta due to excluded solvent pairs...");
  const unsigned int num_cg = topo.num_chargegroups();
  for (unsigned int i = topo.num_solute_chargegroups() + rank; i < num_cg; i += size) {
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
        Storage & storage, int rank, int size) {
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
  for (unsigned int i = 0; i < num_atoms; ++i) {
    periodicity.put_into_positive_box(r(i));

    DEBUG(11, "r(" << i << ") in box: " << math::v2s(r(i)));
  }

  // cache for sin(k . r) and cos(k .  r) terms
  math::SArray sin_kr(num_atoms, 0.0), cos_kr(num_atoms, 0.0);

  // on fly we can calculate the methodology dependent A2 term
  double a2_tilde = 0.0;

  // virial stuff
  const bool do_virial = sim.param().pcouple.virial != math::no_virial;
  math::SymmetricMatrix virial(0.0);
  math::SymmetricMatrix sum_gammahat(0.0);

  // loop over k space
  for (; it != to; ++it) {
    DEBUG(12, "k: " << math::v2s(it->k));
    double C_k = 0.0;
    double S_k = 0.0;
    // loop over atoms, calculate C/S_k and cache the sins and coses
    for (unsigned int i = 0; i < num_atoms; ++i) {
      const double r_dot_k = math::dot(it->k, r(i));
      sin_kr(i) = sin(r_dot_k);
      cos_kr(i) = cos(r_dot_k);
      C_k += topo.charge(i) * cos_kr(i);
      S_k += topo.charge(i) * sin_kr(i);
    }

    // calculate the force component from k
    for (unsigned int i = 0; i < num_atoms; ++i) {
      f(i) += it->k * (it->k2i_gammahat * (C_k * sin_kr(i) - S_k * cos_kr(i)));
    }
    // and the energy component from k
    const double ewald_factor = C_k * C_k + S_k * S_k;
    energy += it->k2i_gammahat * ewald_factor;
    // add to the A2 term
    a2_tilde += it->k2i_gammahat;
    // do the virial
    if (do_virial) {
      const double isotropic_factor =
              (it->ak_gammahat_prime - 2.0 * it->fourier_coefficient) *
              (it->k2i * it->k2i);
      const math::SymmetricMatrix term(math::symmetric_tensor_product(isotropic_factor * it->k, it->k));
      sum_gammahat += term;
      virial += ewald_factor * term;
      virial.add_to_diagonal(ewald_factor * it->k2i_gammahat);
    } // if virial
  } // for k space

  // loop again over atoms and store force
  for (unsigned int i = 0; i < num_atoms; ++i) {
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

  // add the A2 sum before it is scaled
  if (do_virial)
    sum_gammahat.add_to_diagonal(a2_tilde);

  // scale the a2 term
  a2_tilde *= 4.0 * math::Pi / volume;
  conf.lattice_sum().a2_tilde = a2_tilde;

  // and the virial
  if (do_virial) {
    virial *= -0.25 * eps_volume_i;
    DEBUG(6, "Ewald k-space virial:\n" << math::m2s(virial));
// !!! multiply by -2 because of multiplication with -0.5 in pressure calculation
   // storage.virial_tensor += virial;
    storage.virial_tensor += (-2.0 * virial);

    conf.lattice_sum().a2_tilde_derivative = (4.0 * math::Pi / volume) * sum_gammahat;
  }

}

void interaction::Nonbonded_Outerloop
::ls_p3m_kspace_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size,
        util::Algorithm_Timer & timer,
        bool & is_ok) {
  SPLIT_INNERLOOP(_ls_p3m_kspace_outerloop, topo, conf, sim,
          storage, rank, size, timer, is_ok);
}

template<typename t_interaction_spec>
void interaction::Nonbonded_Outerloop
::_ls_p3m_kspace_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage, int rank, int size,
        util::Algorithm_Timer & timer,
        bool & is_ok) {
  DEBUG(7, "\tcalculate interactions in k-space (P3M)");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  math::VArray r = conf.current().pos;
  // Do we have to gather here ?
  const unsigned int num_atoms = topo.num_atoms();
  for (unsigned int i = 0; i < num_atoms; ++i) {
    periodicity.put_into_positive_box(r(i));

    DEBUG(11, "r(" << i << ") in box: " << math::v2s(r(i)));
  }

  if (sim.openmp && rank != 0) {
    // for openmp this is single threaded
    return;
  }
  
  // decompose into domains
  if (sim.mpi)
    interaction::Lattice_Sum::decompose_into_domains<configuration::ParallelMesh > (topo, conf, sim, storage.domain, r, size);
  else
    interaction::Lattice_Sum::decompose_into_domains<configuration::Mesh > (topo, conf, sim, storage.domain, r, size);

  DEBUG(10, "size domain(" << rank << "): " << storage.domain.size());

  // always calculate at the beginning or read it from file.
  if (sim.steps() == 0 && !sim.param().nonbonded.influence_function_read) {
    DEBUG(10, "\t calculating influence function");
    if (rank == 0)
      timer.start("P3M: influence function");
    if (sim.mpi)
      conf.lattice_sum().influence_function.template calculate< configuration::ParallelMesh > (topo, conf, sim);
    else
      conf.lattice_sum().influence_function.template calculate< configuration::Mesh > (topo, conf, sim);

    const double new_rms_force_error = conf.lattice_sum().influence_function.rms_force_error();
    DEBUG(10, "\testimated force error: " << new_rms_force_error);

    if (rank == 0)
      timer.stop("P3M: influence function");

    if (new_rms_force_error > sim.param().nonbonded.influence_function_rms_force_error) {
      io::messages.add("P3M: RMS force error is still too big after reevaluation "
              "of the influence function. Increase the number of grid points, "
              "the rms force error threshold, or the charge width parameter.",
              "P3M", io::message::error);
// !!! this does not lead to a crash of the simulation, but it should crash
// set a bool which is checked in nonbonded_set 
     is_ok = false;
      return;
    }
  }

  // give the current box to the influence function for correction
  // this has to be done AFTER the influence function is calculated
  conf.lattice_sum().influence_function.setBox(conf.current().box);

  // check whether we have to update the influence function
  if (sim.steps() && sim.param().nonbonded.accuracy_evaluation &&
          sim.steps() % sim.param().nonbonded.accuracy_evaluation == 0) {
    if (rank == 0)
      timer.start("P3M: accuracy evaluation");
    if (sim.mpi)
      conf.lattice_sum().influence_function.template evaluate_quality<configuration::ParallelMesh > (topo, conf, sim);
    else
      conf.lattice_sum().influence_function.template evaluate_quality<configuration::Mesh > (topo, conf, sim);
    // number of charges is set to number of atoms as in promd...
    // see MD02.10 eq. C7
    const double rms_force_error = conf.lattice_sum().influence_function.rms_force_error();
    if (rank == 0)
      timer.stop("P3M: accuracy evaluation");

    if (rms_force_error > sim.param().nonbonded.influence_function_rms_force_error) {
      // recalculate the influence function
      DEBUG(10, "\t calculating influence function");
      if (rank == 0)
        timer.start("P3M: influence function");
      if (sim.mpi)
        conf.lattice_sum().influence_function.template calculate<configuration::ParallelMesh > (topo, conf, sim);
      else
        conf.lattice_sum().influence_function.template calculate<configuration::Mesh > (topo, conf, sim);

      const double new_rms_force_error = conf.lattice_sum().influence_function.rms_force_error();
      if (rank == 0)
        timer.stop("P3M: influence function");

      if (new_rms_force_error > sim.param().nonbonded.influence_function_rms_force_error) {
        io::messages.add("P3M: RMS force error is still too big after reevaluation "
                "of the influence function. Increase the number of grid points, "
                "the rms force error threshold, or the charge width parameter.",
                "P3M", io::message::error);
// !!! this does not lead to a crash of the simulation, but it should crash
// set a bool which is checked in nonbonded_set 
     is_ok = false;
        return;
      }
    } // if force error
  } // if update influence function

  // check whether we need to calculate the A2~ term via P3M.
  // check whether we have to calculate it at all in this step
  bool calculate_lattice_sum_corrections =
          sim.param().pcouple.scale != math::pcouple_off || // NPT - every step
          !sim.steps() || // at the beginning of the simulation
          (sim.param().print.stepblock != 0 && sim.steps() % abs(sim.param().print.stepblock) == 0) ||   // energy output req.
          (sim.param().write.energy != 0 && sim.steps() % abs(sim.param().write.energy) == 0); // energy output req. !!!
  const bool do_a2t =
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact_a2_numerical ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_ave_a2_numerical;

  // calculate the A2~ self term via P3M.
  if (do_a2t && calculate_lattice_sum_corrections) {
    DEBUG(10, "\tstarting to assign squared charge to grid ... ");
    if (rank == 0)
      timer.start("P3M: self term");

    if (sim.param().nonbonded.ls_calculate_a2 != simulation::ls_a2t_ave_a2_numerical) {
      // calculate the real A2~ term (not averaged)
      if (sim.mpi)
        interaction::Lattice_Sum::calculate_squared_charge_grid<configuration::ParallelMesh > (topo, conf, sim, storage.domain, r);
      else
        interaction::Lattice_Sum::calculate_squared_charge_grid<configuration::Mesh > (topo, conf, sim, storage.domain, r);
    } else {
      // calculate the averaged A2~ term
      if (sim.mpi)
        interaction::Lattice_Sum::calculate_averaged_squared_charge_grid<configuration::ParallelMesh > (topo, conf, sim, storage.domain);
      else
        interaction::Lattice_Sum::calculate_averaged_squared_charge_grid<configuration::Mesh > (topo, conf, sim, storage.domain);
    }
    DEBUG(10, "\tstarting fft of the squared charge");
    // FFT the charge density grid
    configuration::Mesh & squared_charge = *conf.lattice_sum().squared_charge;
    squared_charge.fft(configuration::Mesh::fft_forward);

    DEBUG(10, "\tcalculation of self term via P3M.");
    if (sim.mpi) {
      interaction::Lattice_Sum::calculate_p3m_selfterm<configuration::ParallelMesh > (topo, conf, sim);
    } else {
      interaction::Lattice_Sum::calculate_p3m_selfterm<configuration::Mesh > (topo, conf, sim);
    }

    if (rank == 0)
      timer.stop("P3M: self term");
  }
  
  
  if (rank == 0)
    timer.start("P3M: energy & force");

  DEBUG(10, "\t done with influence function, starting to assign charge density to grid ... ");
  if (sim.mpi)
    interaction::Lattice_Sum::calculate_charge_density<configuration::ParallelMesh > (topo, conf, sim, storage.domain, r);
  else
    interaction::Lattice_Sum::calculate_charge_density<configuration::Mesh > (topo, conf, sim, storage.domain, r);


  DEBUG(10, "\t assigned charge density to grid, starting fft of charge density");
  // FFT the charge density grid
  configuration::Mesh & charge_density = *conf.lattice_sum().charge_density;
  charge_density.fft(configuration::Mesh::fft_forward);
  DEBUG(10, "\t done with fft! Starting to calculate the energy ...");

  if (sim.mpi) {
    interaction::Lattice_Sum::calculate_potential_and_energy<configuration::ParallelMesh > (topo, conf, sim, storage);
  } else {
    interaction::Lattice_Sum::calculate_potential_and_energy<configuration::Mesh > (topo, conf, sim, storage);
  }
  DEBUG(10, "\t done with calculation of elec. potential and energy, calculating electric field...");
  if (sim.mpi) {
    interaction::Lattice_Sum::calculate_electric_field<configuration::ParallelMesh > (topo, conf, sim);
  } else {
    interaction::Lattice_Sum::calculate_electric_field<configuration::Mesh > (topo, conf, sim);
  }
  DEBUG(10, "\t done with electric field calculation, calculating forces");
  if (sim.mpi) {
    interaction::Lattice_Sum::calculate_force<configuration::ParallelMesh > (topo, conf, sim, storage, r);
  } else {
    interaction::Lattice_Sum::calculate_force<configuration::Mesh > (topo, conf, sim, storage, r);
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
  DEBUG(8, "\telectrostatic self energy and virial");
#ifdef OMP
  if (rank != 0) return;
  size = 1;
#endif

  double a1 = 0.0, a3 = 0.0;
  double & a2_tilde = conf.lattice_sum().a2_tilde;
  double a2 = 0.0;

  const double st2 = topo.sum_squared_charges();
  const double s2 = topo.squared_sum_charges();

  const bool do_virial = sim.param().pcouple.virial != math::no_virial;

  const int shape = sim.param().nonbonded.ls_charge_shape;
  // see MD05.32 Table 6
  switch (shape) {
    case -1:
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

  math::SymmetricMatrix a1_self_term_virial(0.0);
  math::SymmetricMatrix a1_constant_term_virial(0.0);
  // calculate the virial contribution of the A1 term

  if (do_virial) {
    // A1 virial
    const double a1_isotropic_self_virial = -0.25 * math::four_pi_eps_i * a1 * st2;
    DEBUG(10, "\ta1 self virial: " << a1_isotropic_self_virial);
    a1_self_term_virial.add_to_diagonal(a1_isotropic_self_virial);
    const double a1_isotopic_constant_virial = -0.25 * math::four_pi_eps_i * a1 * (s2 - st2);
    DEBUG(10, "\ta1 constant virial: " << a1_isotopic_constant_virial);
    a1_constant_term_virial.add_to_diagonal(a1_isotopic_constant_virial);
  }

  math::SymmetricMatrix a2_virial(0.0);
  // calculate the a2 term
  // do we have to do it numerically?
  if (sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2_numerical ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact_a2_numerical ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_ave_a2_numerical) {

    // check whether the box is a cube.
    bool box_is_cube = sim.param().boundary.boundary == math::rectangular &&
            conf.current().box(0)(0) == conf.current().box(1)(1) &&
            conf.current().box(0)(0) == conf.current().box(2)(2);

    if (box_is_cube) {
      // use precalculated values for A2 and the virial. The cubic case
      DEBUG(10, "a2: using precalculated values for cube.");
      const double wigner_cube = -2.83729748;
      a2 = wigner_cube / conf.current().box(0)(0) - a1 - a3;
      if (do_virial) {
        math::SymmetricMatrix a2_derivative(0.0);
        a2_derivative.add_to_diagonal(wigner_cube / (3.0 * conf.current().box(0)(0)) - a1);
        a2_virial = (-0.25 * st2 * math::four_pi_eps_i) * a2_derivative;
      }
    } else {
      // calculate A2 numerically. This is the rectangular - triclinc case.
      const double & required_precision = sim.param().nonbonded.ls_a2_tolerance;
      DEBUG(10, "\tA2 tolerance: " << required_precision)
      math::Matrix l_to_k = configuration::KSpace_Utils::l_to_k_matrix(
              conf.current().box, conf.boundary_type);

      // Here, we loop over the surface of triclinic volumes of increasing
      // size in l-space, and add the successive A2 contributions of these
      // surfaces.
      math::SymmetricMatrix sum_gammahat(0.0);
      int l_max = 0;
      double tolerance = 0.0;
      a2 = 0.0;
      do {
        ++l_max;
        double term = 0.0;

        // the planes are located perpendicular to the axis (coord)
        for (unsigned int coord = 0; coord < 3; ++coord) {
          unsigned int coord_a = 0, coord_b = 0;
          int boundary_a = 0, boundary_b = 0;
          switch (coord) {
            case 0: // plane perpendicular to x
              coord_a = 1;
              coord_b = 2;
              boundary_a = l_max;
              boundary_b = l_max;
              break;
            case 1: // plane perpendicular to y
              coord_a = 0;
              coord_b = 2;
              // exclude ks in the x plane already considered
              boundary_a = l_max - 1;
              boundary_b = l_max;
              break;
            case 2: // plane perpendicular to z
              coord_a = 0;
              coord_b = 1;
              // exclude ks in the x and y plane already consideres
              boundary_a = l_max - 1;
              boundary_b = l_max - 1;
              break;
          }

          // the plane can be located at -l_max or +l_max but we restrict
          // to summation to the planes localted on +l_max and multiply the resulting
          // term and derivative by a factor of 2.0
          DEBUG(12, "\tnew plane");

          // loop over the plane excluding edges for some axes
          // here we can introduce parallelization by stride
#ifdef OMP
#pragma omp parallel for
#endif
          for (int l_a = -boundary_a + rank; l_a <= boundary_a; l_a += size) {
            math::GenericVec<int> l(0);
            l(coord) = l_max;
            l(coord_a) = l_a;
            for (int l_b = -boundary_b; l_b <= boundary_b; ++l_b) {
              l(coord_b) = l_b;

              DEBUG(13, "\t\tl: " << math::v2s(l));
              const math::Vec & k = math::product(l_to_k, l);
              const double k2 = math::abs2(k);
              const double abs_k = sqrt(k2);
              const double ak = abs_k * width;

              double gamma_hat = 0.0, gamma_hat_prime = 0.0;
              if (do_virial) {
                interaction::Lattice_Sum::charge_shape_fourier(shape,
                        ak, gamma_hat, &gamma_hat_prime);
              } else {
                interaction::Lattice_Sum::charge_shape_fourier(shape,
                        ak, gamma_hat);
              }

#ifdef OMP
#pragma omp critical(addup_term)
#endif
              term += gamma_hat / k2;
              DEBUG(13, "\t\t\tgamma_hat / k2: " << gamma_hat / k2);

              if (do_virial) {
                // factor 2.0 is due to symmetry
                const double isotropic_factor = 2.0 * (ak * gamma_hat_prime - 2.0 * gamma_hat) / (k2 * k2);
#ifdef OMP
#pragma omp critical(addup_sum_gammahat)
#endif
                sum_gammahat += math::symmetric_tensor_product(isotropic_factor * k, k);
              } // virial?
            }
          } // loop over planes         
        } // loop over coordinates

        // now calculate A2 and the relative tolerance
#ifdef XXMPI
        if (sim.mpi) {
          //TODO : CONTORL THAT CORRECT! bschroed
          double my_term = term;
          MPI_Allreduce(&my_term, &term, 1, MPI::DOUBLE, MPI::SUM, sim.mpiControl().comm);
        }
#endif

        // take symmetry into account
        term *= 2.0;

        a2 += term;

        if (do_virial && rank == 0)
          sum_gammahat.add_to_diagonal(term);

        tolerance = fabs(term / a2);
        DEBUG(12, "\ttolerance: " << tolerance);
      } while (tolerance > required_precision);

#ifdef XXMPI
      // for MPI we only have parts of the A2 derivative sum. So we have to add them here
      if (sim.mpi && do_virial) {
        math::SymmetricMatrix sum_gammahat_part = sum_gammahat;

        if (rank) { // slave
          MPI_Reduce(&sum_gammahat_part(0), NULL, 6, MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
        } else { // master
          MPI_Reduce(&sum_gammahat_part(0), &sum_gammahat(0), 6,
                  MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
        }
      }
#endif

      a2 *= 4.0 * math::Pi / volume;
      if (do_virial) {
        const math::SymmetricMatrix a2_derivative = sum_gammahat * (4.0 * math::Pi / volume);
        DEBUG(10, "\tA2 derivative:\n\t" << math::m2s(a2_derivative));
        a2_virial = (-0.25 * st2 * math::four_pi_eps_i) * a2_derivative;
      }
    } // if triclinic

  }
  DEBUG(8, "\tA2 virial:\n" << math::m2s(a2_virial));

#ifdef XXMPI
  // for MPI we only have parts of the A2~ sums. So we have to add them here
  // but only if they were calculated.
  if (sim.mpi && (
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact_a2_numerical ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_ave_a2_numerical)) {

    //TODO: CONTROL THAT CORRECT bscrhoed
    double a2_part = conf.lattice_sum().a2_tilde;
    math::SymmetricMatrix a2_deriv_part = conf.lattice_sum().a2_tilde_derivative;

    if (rank) { // slave    \\TODO: WHY???? bschroed
      MPI_Reduce(&a2_part, NULL, 1, MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&a2_deriv_part(0), NULL, 6, MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
    } else { // master
      MPI_Reduce(&a2_part, &conf.lattice_sum().a2_tilde, 1,
              MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
      MPI_Reduce(&a2_deriv_part(0), &conf.lattice_sum().a2_tilde_derivative(0), 6,
              MPI::DOUBLE, MPI::SUM, sim.mpiControl().masterID, sim.mpiControl().comm);
    }
  }
#endif
  math::SymmetricMatrix a2_tilde_virial(0.0);
  if (do_virial && (
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact_a2_numerical ||
          sim.param().nonbonded.ls_calculate_a2 == simulation::ls_a2t_ave_a2_numerical)) {
    // calculate the virial of the methodology dependent A2 term
    a2_tilde_virial = (-0.25 * st2 * math::four_pi_eps_i) *
            conf.lattice_sum().a2_tilde_derivative;
    DEBUG(8, "A2 tilde virial: " << math::m2s(a2_tilde_virial));
  }

  switch (sim.param().nonbonded.ls_calculate_a2) {
    case simulation::ls_a2_zero:
      // we just set both A2 to zero.
      a2 = a2_tilde = 0.0;
      a2_virial = a2_tilde_virial = 0.0;
      break;
    case simulation::ls_a2t_exact:
      // A2t was calculated exactly by Ewald/P3M. A2 is set to A2t
      a2 = a2_tilde;
      a2_virial = a2_tilde_virial;
    case simulation::ls_a2_numerical:
      // A2 was calculate numerically, A2t is set to A2
      a2_tilde = a2;
      a2_tilde_virial = a2_virial;
      break;
    case simulation::ls_a2t_exact_a2_numerical:
    case simulation::ls_a2t_ave_a2_numerical:
      // we already have A2t and A2 and the virials
      break;
    default:
      io::messages.add("A2 calculation method not implemented", "Lattice Sum",
              io::message::critical);
  } // switch ls_calculate_a2

  DEBUG(8, "a1 = " << a1 << ", a2 = " << a2 << ", ~a2 = " << a2_tilde << ", a3 = " << a3);
  // now combine everything (MD05.32 eq 54)
  if (rank == 0) {
    storage.energies.ls_self_total = (a1 + a2 + a3) * st2 * math::eps0_i / (8.0 * math::Pi);
    if (do_virial) {
      // we have to remove the A2 virial from the self term virial
      const math::SymmetricMatrix self_term_virial = a1_self_term_virial + a2_virial;
// !!! since multiplication with -0.5 in pressure calc, multiply with -2 here
     // storage.virial_tensor += self_term_virial;
      storage.virial_tensor +=  ( -2.0 * self_term_virial );
      DEBUG(6, "\tself term virial:\n\t" << math::m2s(self_term_virial));
    } // virial
  } else {
    storage.energies.ls_self_total = 0.0;
  }

  DEBUG(6, "ls_self_total = " << storage.energies.ls_self_total);

  // now calculate the E_A term
  if (rank == 0) {
    storage.energies.ls_a_term_total = (a1 * s2 - (a1 + a2_tilde) * st2) * math::eps0_i / (8.0 * math::Pi);
    if (do_virial) {
      // We have to add the A2~ virial to the constant term virial
      const math::SymmetricMatrix constant_term_virial = a1_constant_term_virial - a2_tilde_virial;
// !!! since multiplication with -0.5 in pressure calc, multiply with -2 here
     // storage.virial_tensor += constant_term_virial;
      storage.virial_tensor += ( -2.0 *  constant_term_virial );
      DEBUG(6, "\tconstant term virial:\n\t" << math::m2s(constant_term_virial));
    } // virial
  } else {
    storage.energies.ls_a_term_total = 0.0;
  }
  DEBUG(6, "ls_a_term_total = " << storage.energies.ls_a_term_total);

  if (sim.param().pcouple.scale == math::pcouple_off && !sim.steps()){ // NVT and first step
    conf.current().energies.ls_self_total_nvt = storage.energies.ls_self_total;
    conf.old().energies.ls_self_total_nvt = storage.energies.ls_self_total;
    conf.current().energies.ls_a_term_total_nvt = storage.energies.ls_a_term_total;
    conf.old().energies.ls_a_term_total_nvt = storage.energies.ls_a_term_total;
  }
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
        Storage & storage, int rank, int size) {
  DEBUG(8, "\telectrostatic surface terms");

  // under tinfoil boundary conditions we don't have
  // to calculate anything as the terms are zero.
  if (sim.param().nonbonded.ls_epsilon == 0.0) {
    DEBUG(10, "\ttinfoil.");
    storage.energies.ls_surface_total = 0.0;
    return;
  }

  math::Vec box_dipole_moment(0.0);
  math::Vec box_centre = conf.current().box(0) / 2.0 +
          conf.current().box(1) / 2.0 +
          conf.current().box(2) / 2.0;

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);

  const unsigned int num_atoms = topo.num_atoms();
  for (unsigned int i = 0; i < num_atoms; i++) {
    math::Vec r = conf.current().pos(i);
    box_dipole_moment += topo.charge(i) * (r - box_centre);
  }

  // now we have the box dipole moment and can calculate the energy and force
  const double prefactor = math::eps0_i /
          ((sim.param().nonbonded.ls_epsilon * 2.0 + 1.0) *
          math::volume(conf.current().box, conf.boundary_type));

  const double abs2_box_dipole_moment = math::abs2(box_dipole_moment);
  storage.energies.ls_surface_total = 0.5 * abs2_box_dipole_moment * prefactor;
  DEBUG(6, "\tsurface energy: " << storage.energies.ls_surface_total);

  for (unsigned int i = 0; i < num_atoms; i++) {
    storage.force(i) += -prefactor * topo.charge(i) * box_dipole_moment;
  }

  // do the virial
  if (sim.param().pcouple.virial != math::no_virial) {
    math::Matrix virial(0.0);
    const double isotropic_term = -0.25 * prefactor * abs2_box_dipole_moment;
    for (unsigned int i = 0; i < 3; ++i) {
      virial(i, i) = isotropic_term
              + 0.5 * prefactor * box_dipole_moment(i) * box_dipole_moment(i);
    }
    storage.virial_tensor += virial;
    DEBUG(6, "surface term virial: " << math::m2s(virial));
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
        ) {
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
        double & e_lj, double & e_crf) {
  math::Vec r;
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);

  Nonbonded_Term term;
  term.init(sim);

  const lj_parameter_struct &lj =
          m_param.lj_parameter(topo.iac(atom_i),
          topo.iac(atom_j));

  periodicity.nearest_image(conf.current().pos(atom_i), conf.current().pos(atom_j), r);
  DEBUG(10, "\tni i " << conf.current().pos(atom_i)(0) << " / "
          << conf.current().pos(atom_i)(1) << " / "
          << conf.current().pos(atom_i)(2));
  DEBUG(10, "\tni j " << conf.current().pos(atom_j)(0) << " / "
          << conf.current().pos(atom_j)(1) << " / "
          << conf.current().pos(atom_j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  double f = 0.0;
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
        PairlistContainer const & pairlist) {

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
        PairlistContainer const & pairlist) {

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

    for (unsigned int d1 = 0; d1 < 3; ++d1)
      for (unsigned int d2 = 0; d2 < 3; ++d2)
        hessian(d1, d2) += h(d1, d2);
  }

  return 0;
}
