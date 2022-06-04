/**
 * @file berendsen_barostat.cc
 * methods of the berendsen barostat.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../configuration/state_properties.h"

#include "../../util/error.h"

#include "berendsen_barostat.h"
#include "../../math/boundary_checks.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE pressure

#include "../../util/debug.h"
#include "../../math/transformation.h"

int algorithm::Berendsen_Barostat
::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  m_timer.start();

  DEBUG(8, "Berendsen Barostat == apply");

  // position are current!
  math::VArray & pos = conf.current().pos;
  math::VArray & ref = conf.special().reference_positions;
  math::Matrix & pressure = conf.old().pressure_tensor;
  math::Box & box = conf.current().box;

  DEBUG(9, "scaling = " << sim.param().pcouple.scale);

  bool scale_ref = sim.param().posrest.posrest != simulation::posrest_off &&
          sim.param().posrest.scale_reference_positions;

  switch (sim.param().pcouple.scale) {
    case math::pcouple_isotropic :
    {
      DEBUG(9, "total pressure (isotropic) ...");
      double total_pressure = (pressure(0, 0)
              + pressure(1, 1)
              + pressure(2, 2)) / 3.0;

      DEBUG(8, "pressure: " << total_pressure);

      double mu = pow(1.0 - sim.param().pcouple.compressibility
              * sim.time_step_size() / sim.param().pcouple.tau
              * (sim.param().pcouple.pres0(0, 0) - total_pressure),
              1.0 / 3.0);

      DEBUG(8, "mu: " << mu);

      // scale the box
      box *= mu;

      // scale the positions
      for (unsigned int i = 0; i < pos.size(); ++i) {
        pos(i) = mu * pos(i);
        if (scale_ref)
          ref(i) = mu * ref(i);
      }
      break;
    }
    case math::pcouple_anisotropic :
    {
      math::Vec mu;

      DEBUG(8, "anisotropic pressure scaling");

      for (int i = 0; i < 3; ++i) {
        mu(i) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * (sim.param().pcouple.pres0(i, i) -
                pressure(i, i)),
                1.0 / 3.0);
      }

      DEBUG(10, "mu = " << math::v2s(mu));

      // scale the box
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          box(i)(j) *= mu(j);

      DEBUG(10, "and the positions...");

      // scale the positions
      for (unsigned int i = 0; i < pos.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
          pos(i)(j) *= mu(j);
          if (scale_ref)
            ref(i)(j) *= mu(j);
        }
      }
      break;
    }
    case math::pcouple_full_anisotropic :
    {

      math::Matrix mu;
      double delta = 0.0, mu_aux = 0.0;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          if (i == j)
            delta = 1;
          else
            delta = 0;
          mu_aux = delta - sim.param().pcouple.compressibility
                  * sim.time_step_size() / (3*sim.param().pcouple.tau)
                  * (sim.param().pcouple.pres0(i, j) -
                  pressure(i, j));
          mu(i, j) = mu_aux;
     //     mu(i, j) = math::sign(mu_aux) * pow(fabs(mu_aux),
      //            1.0 / 3.0);

        }
      }

      // scale the box
      // decompose in rotation and boxlength/boxangle scaling
    /*  math::Matrix Rmu(math::rmat(mu));
      math::Matrix mu_new = math::product(math::transpose(Rmu), mu);
      DEBUG(10, "mu: \n" << math::m2s(mu));
      DEBUG(10, "Rmu: \n" << math::m2s(Rmu));
      DEBUG(10, "mu_new: \n" << math::m2s(mu_new));*/
      
      
      //Bruno & Pitschna: symmetrise as in GROMACS (29 Sep, 2010)
      for (int i = 0; i < 3; ++i) 
        for (int j = 0; j < 3; ++j) 
          if (i<j)
            mu(i, j) = mu(i, j) + mu(j, i);
        
      
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          if (i>j)
            mu(i, j) = 0;
          //std::cout << "  " << mu(i, j) << "  ";
          
        }
        //std::cout << "\n";
      }

      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          //std::cout << "  " << box(i, j) << "  ";
        }
        //std::cout << "\n";
      }


        
      
      //Scale the box
      math::Box m(0.0);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 3; ++k)
            m(i)(j) += mu(k, i) * box(k)(j);
            //m(j)(i) += mu(i, k) * box(j)(k);
      box = m;

   
      //we want a triangular matrix
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if (fabs(box(j)(i)) < math::epsilon)
            box(j)(i) = 0;

      DEBUG(10, "new box: \n" << math::m2s(math::Matrix(box)));
      //scale the positions
      for (unsigned int i = 0; i < pos.size(); ++i) {
        pos(i) = math::product(mu, pos(i));
        if (scale_ref) {
          ref(i) = math::product(mu, ref(i));
        }
      }


      // new Euleur angles
      /* neglect the rotational contribution for now... (gives weird results at the moment :o)...)
       * DEBUG(10, "old Euler angles: " << conf.current().phi
                 << " " << conf.current().theta
                 << " " << conf.current().psi);
        math::Matrixl Rmat(math::rmat(conf.current().phi, conf.current().theta, conf.current().psi));
        math::Matrixl RmatRmu(math::product(Rmat, Rmu));
        long double R11R21 = sqrtl(RmatRmu(0, 0) * RmatRmu(0, 0) + RmatRmu(0, 1) * RmatRmu(0, 1));
        if (R11R21 == 0.0) {
          conf.current().theta = -math::sign(RmatRmu(0, 2)) * M_PI / 2;
          conf.current().psi = 0.0;
          conf.current().phi = -math::sign(RmatRmu(1, 0)) * acosl(math::costest(RmatRmu(1, 1)));
        } else {
          conf.current().theta = -math::sign(RmatRmu(0, 2)) * acosl(math::costest(R11R21));
          long double costheta = cosl(conf.current().theta);
          conf.current().psi = math::sign(RmatRmu(1, 2) / costheta) * acosl(math::costest(RmatRmu(2, 2) / costheta));
          conf.current().phi = math::sign(RmatRmu(0, 1) / costheta) * acosl(math::costest(RmatRmu(0, 0) / costheta));
        }
         conf.old().phi=conf.current().phi;
         DEBUG(10, "new Euler angles: " << conf.current().phi
                 << " " << conf.current().theta
                 << " " << conf.current().psi);
       */
    }
    case math::pcouple_semi_anisotropic :
    {
      math::Vec mu;

      DEBUG(8, "semi anisotropic pressure scaling");

      if (sim.param().pcouple.x_semi < 1)
        mu(0) = 1;
      else{
        mu(0) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * (sim.param().pcouple.pres0(0, 0) -
                pressure(0, 0)),
                1.0 / 3.0);
      }

      if (sim.param().pcouple.y_semi < 1)
        mu(1) = 1;
      if (sim.param().pcouple.y_semi > 0 && (sim.param().pcouple.y_semi == sim.param().pcouple.x_semi)){
        mu(0) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * ((sim.param().pcouple.pres0(0, 0) + sim.param().pcouple.pres0(1, 1) -
                pressure(0, 0) - pressure(1, 1))/2),
                1.0 / 3.0);
        mu(1) = mu(0);
      }
      if (sim.param().pcouple.y_semi > 0 && (sim.param().pcouple.y_semi != sim.param().pcouple.x_semi)){
        mu(1) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * (sim.param().pcouple.pres0(1, 1) -
                pressure(1, 1)),
                1.0 / 3.0);
      }

      if(sim.param().pcouple.z_semi < 1)
        mu(2) = 1;
      if(sim.param().pcouple.z_semi > 0) {
        if(sim.param().pcouple.z_semi == sim.param().pcouple.x_semi){
          mu(0) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * ((sim.param().pcouple.pres0(0, 0) + sim.param().pcouple.pres0(2, 2) -
                pressure(0, 0) - pressure(2, 2))/2),
                1.0 / 3.0);
          mu(2) = mu(0);
        }
        if(sim.param().pcouple.z_semi == sim.param().pcouple.y_semi){
          mu(1) = pow(1.0 - sim.param().pcouple.compressibility
                * sim.time_step_size() / sim.param().pcouple.tau
                * ((sim.param().pcouple.pres0(1, 1) + sim.param().pcouple.pres0(2, 2) -
                pressure(1, 1) - pressure(2, 2))/2),
                1.0 / 3.0);
          mu(2) = mu(1);
        }
        if((sim.param().pcouple.z_semi != sim.param().pcouple.x_semi) && (sim.param().pcouple.z_semi != sim.param().pcouple.y_semi)) {
          mu(2) = pow(1.0 - sim.param().pcouple.compressibility
                  * sim.time_step_size() / sim.param().pcouple.tau
                  * (sim.param().pcouple.pres0(2, 2) -
                  pressure(2, 2)),
                  1.0 / 3.0);
        }
      }

      

      DEBUG(10, "mu = " << math::v2s(mu));

      // scale the box
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          box(i)(j) *= mu(j);

      DEBUG(10, "and the positions...");

      // scale the positions
      for (unsigned int i = 0; i < pos.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
          pos(i)(j) *= mu(j);
          if (scale_ref)
            ref(i)(j) *= mu(j);
        }
      }
      break;
    }
    default:
      return 0;
  }

  if (!sim.param().multicell.multicell && !math::boundary_check_cutoff(conf.current().box,
      sim.param().boundary.boundary, sim.param().pairlist.cutoff_long)) {
    io::messages.add("box is too small: not twice the cutoff!",
            "Berendsen_Barostat", io::message::error);
  }

  m_timer.stop();

  return 0;

}
