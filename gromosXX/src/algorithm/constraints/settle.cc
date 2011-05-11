/**
 * @file settle.cc
 * contains the template methods for
 * the class Settle.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/shake.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"

#include "settle.h"
#include <vector>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

int algorithm::Settle::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {

  if (!quiet) {
    os << "SETTLE\n"
       << "\tsolvent\n";
  }

  // check for 3 site water model like SPC, TIP3
  
  // we need a solvent
  if (topo.num_solvents() != 1) {
    io::messages.add("SETTLE does only work if 1 solvent.",
            "Settle", io::message::error);
    return 1;
  }
  
  // we need 3 atoms
  if (topo.solvent(0).num_atoms() != 3) {
    io::messages.add("SETTLE does only work with water like molecules (3 atoms).",
            "Settle", io::message::error);    
    return 1;
  }
  
  // the masses of the second and third atom have to be the same
  if (topo.solvent(0).atom(1).mass != topo.solvent(0).atom(2).mass) {
    io::messages.add("SETTLE does only work with water like molecules (wrong masses).",
            "Settle", io::message::error);    
    return 1;
  }

  // the molecule must be rigid: 3 distance constraints
  if (topo.solvent(0).distance_constraints().size() != 3) {
    io::messages.add("SETTLE does only work with water like molecules (3 distance constraints).",
            "Settle", io::message::error);
    return 1;
  }
  
  // the molecule must have two equal bond lengths (constraints 1 and 2)
  if (parameter()[topo.solvent(0).distance_constraint(0).type].r0 !=
      parameter()[topo.solvent(0).distance_constraint(1).type].r0) {
    io::messages.add("SETTLE does only work with water like molecules (distance constraints wrong).",
            "Settle", io::message::error);
    return 1;
  }
  
  // insert the solvent atoms to the constrained atoms
  for(unsigned int i = topo.num_solvent_atoms(); i < topo.num_atoms(); ++i) {
    constrained_atoms().insert(i);
  }
  
  // check whether we do an initial apply of the constraint algorithm
  if (sim.param().start.shake_pos || sim.param().start.shake_vel) {
    io::messages.add("initial settle-ing is not possible.",
            "Settle", io::message::error);
  }
  
  if (!quiet) {
    os << "END\n";
  }
  
  return 0;
}

int algorithm::Settle::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(7, "applying SETTLE");
  m_timer.start();

  int error = 0;

  if (sim.param().system.nsm &&
      sim.param().constraint.solvent.algorithm == simulation::constr_settle) {

    DEBUG(8, "\twe need to settle SOLVENT");
    solvent(topo, conf, sim, error);
    
    if (error) {
      std::cout << "SETTLE: exiting with error condition "
              << "at step " << sim.steps() << std::endl;
      io::messages.add("SETTLE error", "Settle", io::message::error);
      conf.special().shake_failure_occurred = true;
      m_timer.stop();
      return 1;
    }
  }

  m_timer.stop();
  // return success!
  return 0;
}

/** 
 * a helper function to raise a settle error
 */
void algorithm::Settle
::printError(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int atom, std::string message) {
  unsigned int mol = (atom - topo.num_solute_atoms()) / topo.solvent(0).num_atoms() + 1;
  std::cout << "SETTLE ERROR\n"
          << "\tfailed to settle solvent molecule " << mol << ":\n"
          << "\treason: " << message << "\n"
          << "\tref OW:     " << math::v2s(conf.old().pos(atom)) << "\n"
          << "\tref H1:     " << math::v2s(conf.old().pos(atom + 1)) << "\n"
          << "\tref H2:     " << math::v2s(conf.old().pos(atom + 2)) << "\n"
          << "\tpos OW:     " << math::v2s(conf.current().pos(atom)) << "\n"
          << "\tpos H1:     " << math::v2s(conf.current().pos(atom + 1)) << "\n"
          << "\tpos H2:     " << math::v2s(conf.current().pos(atom + 2)) << "\n"
          << "\tvel OW:     " << math::v2s(conf.current().vel(atom)) << "\n"
          << "\tvel H1:     " << math::v2s(conf.current().vel(atom + 1)) << "\n"
          << "\tvel H2:     " << math::v2s(conf.current().vel(atom + 2)) << "\n"
          << "\told vel OW: " << math::v2s(conf.old().vel(atom)) << "\n"
          << "\told vel H1: " << math::v2s(conf.old().vel(atom + 1)) << "\n"
          << "\told vel H2: " << math::v2s(conf.old().vel(atom + 2)) << "\n"
          << "\tforce OW:   " << math::v2s(conf.old().force(atom)) << "\n"
          << "\tforce H1:   " << math::v2s(conf.old().force(atom + 1)) << "\n"
          << "\tforce H2:   " << math::v2s(conf.old().force(atom + 2)) << "\n";
}


void algorithm::Settle
::solvent(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim, int & error) {
  DEBUG(8, "\tSETTLE-ing SOLVENT");

  // assume everything is fine
  error = 0;

  assert(topo.num_solvents() == 1);
  assert(topo.solvent(0).num_atoms() == 3);
  // masses
  const double mass_O = topo.solvent(0).atom(0).mass;
  assert(topo.solvent(0).atom(1).mass == topo.solvent(0).atom(2).mass);
  const double mass_H = topo.solvent(0).atom(1).mass;
  // distance constraints
  assert(topo.solvent(0).distance_constraints().size() == 3);
  assert(parameter()[topo.solvent(0).distance_constraints()[0].type].r0 ==
      parameter()[topo.solvent(0).distance_constraints()[1].type].r0);
  const double dist_OH = parameter()[topo.solvent(0).distance_constraints()[0].type].r0;
  const double dist_HH = parameter()[topo.solvent(0).distance_constraints()[2].type].r0;

  // calculate the coordinates of the canonical triangle
  // see Figure 2 (a)
  const double half_mO_div_mH = 0.5 * mass_O / mass_H;
  const double rc = 0.5 * dist_HH;
  const double ra = sqrt(dist_OH * dist_OH - rc * rc) / (1.0 + half_mO_div_mH);
  const double rb = half_mO_div_mH * ra;

  // needed for constraint force, velocity correction and virial
  const double dt_i = 1.0 / sim.time_step_size();
  const double dt2_i = dt_i * dt_i;
  
  const bool do_velocity = !sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze;

  const int num_atoms = topo.num_atoms();

  // loop over all SPC molecules
#ifdef OMP
#pragma omp parallel for
#endif
  for (int i = topo.num_solute_atoms(); i < num_atoms; i += 3) {
    DEBUG(3, "Molecule: " << (i - topo.num_solute_atoms())/3);
    // get the new positions
    math::Vec * pos_new = &conf.current().pos(i);
    // get the old positions
    const math::Vec * const pos_old = &conf.old().pos(i);
    // get the constraint force
    math::Vec * cons_force = &conf.old().constraint_force(i);

    // See the Figures 1 and 2 for a picture.
    
    // vectors in the plane of the old positions
    math::Vec b0 = pos_old[1] - pos_old[0];
    math::Vec c0 = pos_old[2] - pos_old[0];

    // centre of mass of new positions
    const math::Vec d0 = (pos_new[0] * mass_O + 
                          pos_new[1] * mass_H + pos_new[2] * mass_H) /
                         (mass_O + mass_H + mass_H);

    DEBUG(3, "COM= " << math::v2s(d0));

    // move the origin to the centre of mass
    math::Vec a1 = pos_new[0] - d0;
    math::Vec b1 = pos_new[1] - d0;
    math::Vec c1 = pos_new[2] - d0;

    DEBUG(3, "a1= " << math::v2s(a1) << " b1= " << math::v2s(b1) << " c1= " << math::v2s(c1));

    // vectors describing transformation from original coordinate system to
    // the centre of mass originated coordinate system
    math::Vec n0 = math::cross(b0, c0);
    math::Vec n1 = math::cross(a1, n0);
    math::Vec n2 = math::cross(n0, n1);
    n0 = n0 / math::abs(n0); // this can give a NaN but it is very unlikely.
    n1 = n1 / math::abs(n1);
    n2 = n2 / math::abs(n2);

    DEBUG(3, "n0= " << math::v2s(n0) << " n1= " << math::v2s(n1) << " n2= " << math::v2s(n2));

    // generate new normal vectors from the transpose in order to undo
    // the transformation
    const math::Vec m1(n1(0), n2(0), n0(0));
    const math::Vec m2(n1(1), n2(1), n0(1));
    const math::Vec m0(n1(2), n2(2), n0(2));

    // do the transformation to the centre of mass originated coordinate system
    // of the old positions
    b0 = math::Vec(math::dot(n1, b0), math::dot(n2, b0), math::dot(n0, b0));
    c0 = math::Vec(math::dot(n1, c0), math::dot(n2, c0), math::dot(n0, c0));

    // and of the new positions
    a1 = math::Vec(math::dot(n1, a1), math::dot(n2, a1), math::dot(n0, a1));
    b1 = math::Vec(math::dot(n1, b1), math::dot(n2, b1), math::dot(n0, b1));
    c1 = math::Vec(math::dot(n1, c1), math::dot(n2, c1), math::dot(n0, c1));

    DEBUG(3, "rotated: a1= " << math::v2s(a1) << " b1= " << math::v2s(b1) << " c1= " << math::v2s(c1));

    // now we can compute positions of canonical water 
    // the cos is calculate from the square of the sin via sin^2 + cos^2 = 1
    const double sinphi = a1(2) / ra; // this is (A8)
    const double one_minus_sinphi2 = 1.0 - sinphi*sinphi;
    if (one_minus_sinphi2 < 0.0) {
      printError(topo, conf, sim, i, "sin(phi) > 1.0");
      error = 1;
      continue;
    }
    const double cosphi = sqrt(one_minus_sinphi2);

    const double sinpsi = (b1(2) - c1(2)) / (2.0 * rc * cosphi); // this is (A9)
    const double one_minus_sinpsi2 = 1.0 - sinpsi*sinpsi;
    if (one_minus_sinpsi2 < 0.0) {
      printError(topo, conf, sim, i, "sin(psi) > 1.0");
      error = 1;
      continue;
    }
    const double cospsi = sqrt(one_minus_sinpsi2);

    // these are just to avoid recalculations
    const double minus_rb_cosphi = -rb * cosphi;
    const double rc_cospsi = rc * cospsi;
    const double rc_sinpsi_sinphi = rc * sinpsi*sinphi;
    const double rc_sinpsi_cosphi = rc * sinpsi*cosphi;
    
    // do the X. this is (A3)
    const double x_a2 = 0.0;
    const double x_b2 = - rc_cospsi;
    const double x_c2 = rc_cospsi;
    
    // do the Y. this is (A3)
    // I think there is an error in the paper. ra was confused with rc
    const double y_a2 = ra * cosphi; 
    const double y_b2 = minus_rb_cosphi - rc_sinpsi_sinphi;
    const double y_c2 = minus_rb_cosphi + rc_sinpsi_sinphi;
    
    // do the Z components
    const double z_a2 = ra * sinphi; // this is (A5)
    const double z_b2 = -rb * sinphi + rc_sinpsi_cosphi; // this is (A6)
    const double z_c2 = -rb * sinphi - rc_sinpsi_cosphi; // this is (A7)
    
    // now construct the a2, b2 and c2 vectors
    const math::Vec a2(x_a2, y_a2, z_a2);
    const math::Vec b2(x_b2, y_b2, z_b2);
    const math::Vec c2(x_c2, y_c2, z_c2);

    DEBUG(3, "a2= " << math::v2s(a2) << " b2= " << math::v2s(b2) << " c2= " << math::v2s(c2));

    // there are no a0 terms because we've already subtracted the term off 
    // when we first defined b0 and c0.
    // this is (A15) but the equation is not really numbered...
    const double alpha = b2(0) * (b0(0) - c0(0)) + b0(1) * b2(1) + c0(1) * c2(1);
    const double beta = b2(0) * (c0(1) - b0(1)) + b0(0) * b2(1) + c0(0) * c2(1);
    const double gamma = b0(0) * b1(1) - b1(0) * b0(1) + c0(0) * c1(1) - c1(0) * c0(1);

    const double alpha2_beta2 = alpha * alpha + beta * beta;
    // this is (A17)
    const double sintheta = (alpha * gamma -
            beta * sqrt(alpha2_beta2 - gamma * gamma)) / alpha2_beta2; 
    const double one_minus_sintheta2 = 1.0 - sintheta * sintheta;
    if (one_minus_sintheta2 < 0.0) {
      printError(topo, conf, sim, i, "sin(theta) > 1.0");
      error = 1;
      continue;
    }
    const double costheta = sqrt(one_minus_sintheta2);

    // new finally a3, b3 and c3. These are (A4)
    const math::Vec a3(-a2(1) * sintheta,
            a2(1) * costheta,
            a1(2));
    const math::Vec b3(b2(0) * costheta - b2(1) * sintheta,
            b2(0) * sintheta + b2(1) * costheta,
            b1(2));
    const math::Vec c3(-b2(0) * costheta - c2(1) * sintheta,
            -b2(0) * sintheta + c2(1) * costheta,
            c1(2));

    DEBUG(3, "a3= " << math::v2s(a3) << " b3= " << math::v2s(b3) << " c3= " << math::v2s(c3));
    
    // calculate the new positions
    const math::Vec pos_a = math::Vec(math::dot(a3, m1), math::dot(a3, m2), math::dot(a3, m0)) + d0;
    const math::Vec pos_b = math::Vec(math::dot(b3, m1), math::dot(b3, m2), math::dot(b3, m0)) + d0;
    const math::Vec pos_c = math::Vec(math::dot(c3, m1), math::dot(c3, m2), math::dot(c3, m0)) + d0;

    DEBUG(3, "pos_a= " << math::v2s(pos_a) << " pos_b= " << math::v2s(pos_b) << " pos_c= " << math::v2s(pos_c));

    // calculate the displacement for velocity and virial calculation
    const math::Vec d_a = pos_a - pos_new[0];
    const math::Vec d_b = pos_b - pos_new[1];
    const math::Vec d_c = pos_c - pos_new[2];
    
    // undo the transformation for the a3, b3 and c3 vectors in order to get 
    // the new positions
    pos_new[0] = pos_a;
    pos_new[1] = pos_b;
    pos_new[2] = pos_c;
    
    // in any case calculate the constraint force - it is a very interetsing
    // quantity 
    // by finite difference
    cons_force[0] = d_a * mass_O * dt2_i;
    cons_force[1] = d_b * mass_H * dt2_i;
    cons_force[2] = d_c * mass_H * dt2_i;

    DEBUG(3, "constrained forces 0= " << math::v2s(cons_force[0]) << " 1= " << math::v2s(cons_force[1]) << " 2= " << math::v2s(cons_force[2]));
    
    if (do_velocity) {
      math::Vec * vel_new = &conf.current().vel(i);
      
      // let's reset the velocity
      // recalculate velocity from displacement and timestep
      vel_new[0] += d_a * dt_i;
      vel_new[1] += d_b * dt_i;
      vel_new[2] += d_c * dt_i;
    } // do velocity

    if (sim.param().pcouple.virial == math::atomic_virial) {
      // calculate the constraint virial
#ifdef OMP
#pragma omp critical
#endif
      {
        for (int a = 0; a < 3; ++a) {
          for (int aa = 0; aa < 3; ++aa) {
            conf.old().virial_tensor(a, aa) +=
                    pos_old[0](a) * cons_force[0](aa) +
                    pos_old[1](a) * cons_force[1](aa) +
                    pos_old[2](a) * cons_force[2](aa);
          }
        }
      } // do virial
    }
  } // for molecules

  return;
}
