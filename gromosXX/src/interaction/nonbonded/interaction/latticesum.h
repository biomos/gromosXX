/**
 * @file latticesum.h
 * contains functions for Ewald and P3M
 */

#ifndef INCLUDED_EWALD_H
#define INCLUDED_EWALD_H

namespace interaction {
  class Storage;

  /**
   * @class Lattice_Sum
   * holds static functions for Ewald and P3M. 
   * this class is abstract and you should not create instances of it.
   */
  class Lattice_Sum {
  public:     
    /**
     * the virtual function to ensure the class is not instantiated.
     */
    virtual void dummy() = 0;
    
    /**
     * calculates charge shaping switch function @f$ \eta(\xi) @f$ and its derivative
     * @f$ \eta'(\xi) @f$ 
     * @param[in] shape the charge shape switch
     * @param[out] eta_xi the shaping switch function @f$ \eta(\xi) @f$
     * @param[out] eta_xi_d_xi the derivative @f$ \eta'(\xi) @f$ 
     */
    static inline void charge_shape_switch(const int &shape,
            const double &xi,
            double &eta_xi,
            double &d_eta_xi_d_xi) {
      
      // to do: optimize these functions for speed
      switch(shape) {
        case -1 :
          eta_xi = erfc(xi);
          // - 2*exp(- xi^2) / sqrt(pi) 
          d_eta_xi_d_xi = exp(-xi * xi) * -2.0 / sqrt(math::Pi);
          break;
        case 0 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = 0.5 * pow(1.0 - xi,2.0)*(xi + 2.0) * h1mxi;
          d_eta_xi_d_xi = 3.0 / 2.0 * (xi * xi - 1.0) * h1mxi;
          break;
        }
        case 1 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = pow(1.0 - xi, 3.0) * (xi + 1) * h1mxi;
          d_eta_xi_d_xi = -2.0 * pow(xi - 1.0, 2.0) * (1.0 + 2.0 * xi) * h1mxi;
          break;
        } 
        case 2 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = 1.0 / 2.0 * pow(1.0 - xi, 4.0) * (3.0 * xi + 2.0) * h1mxi;
          d_eta_xi_d_xi = 5.0 / 2.0 * pow(xi - 1.0, 3.0) * (1.0 + 3.0 * xi) * h1mxi;
          break;
        }
        case 3 : {
          const double xi2= xi*xi;
          const double xim1 = xi - 1;
          const double h1mxi = math::heaviside(1-xi);
          eta_xi = 0.25*pow(xim1,4.0)*(4*xi2+7*xi+4)*h1mxi;
          d_eta_xi_d_xi = pow(xim1,3.0)*(2.25+6.75*xi+6*xi2)*h1mxi;
          break;
        }
        case 4 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = 1.0 / 8.0 * pow(1.0 - xi, 5.0) * (15.0 * xi * xi + 19.0 * xi + 8.0) * h1mxi;
          d_eta_xi_d_xi = - 21.0 / 8.0 * pow(xi - 1.0, 4.0) * (1.0 + xi * ( 4.0 + 5.0 * xi)) * h1mxi;
          break;
        }
        case 5 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = pow(1.0 - xi, 6.0) * (3.0 * xi * xi + 3.0 * xi + 1.0) * h1mxi;
          d_eta_xi_d_xi = 3.0 * pow(xi - 1.0, 5.0) * (1.0 + xi * ( 5.0 + 8.0 * xi)) * h1mxi;
          break;
        }
        case 6 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = 1.0 / 16.0 * pow(1.0 - xi, 6.0) * (35.0 * pow(xi,3.0) + 66.0 * xi * xi + 51.0 * xi + 16.0) * h1mxi;
          d_eta_xi_d_xi = 9.0 / 16.0 * pow(xi - 1.0, 5.0) * (5.0 + xi * ( 25.0 + xi * (47.0 + 35.0 * xi ))) * h1mxi;
          break;
        }
        case 7 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = 1.0 / 8.0 * pow(1.0 - xi , 7.0) * (32.0 * pow(xi, 3.0) + 49.0 * xi * xi + 31.0 * xi + 8.0) * h1mxi;
          d_eta_xi_d_xi = -5.0 / 8.0 * pow(xi - 1.0, 6.0) * (5.0 + xi * (30.0 + xi * (69.0 + 64.0 * xi))) * h1mxi;
          break;
        }
        case 8 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = 1.0 / 16.0 * pow(1.0 - xi, 8.0) * (105.0 * pow(xi,3.0) + 136.0 * xi * xi + 73.0 * xi + 16.0) * h1mxi;
          d_eta_xi_d_xi = 55.0 / 16.0 * pow(xi - 1.0, 7.0) * (1.0 + 3.0 * xi) * (1.0 + xi  * (4.0 + 7.0 * xi )) * h1mxi;
          break;
        }
        case 9 : {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi = 1.0 / 32.0 * pow(1.0 - xi, 8.0) * (160.0 * pow(xi, 4.0) + 335.0 * pow(xi,3.0) + 312.0 * xi * xi + 151.0 * xi + 32.0) * h1mxi;
          d_eta_xi_d_xi = 15.0 / 32.0 * pow(xi - 1.0, 7.0) * (7.0 + xi * (49.0 + xi * (141.0 + xi * (203.0 + 128.0*xi)))) * h1mxi;
          break;
        }
        case 10: {
          const double h1mxi = math::heaviside(1.0-xi);
          eta_xi =  1.0 / 128.0 * pow(1.0 - xi, 9.0) * (1155.0 * pow(xi,4.0) + 2075.0 * pow(xi,3.0) + 1665.0 * xi * xi + 697.0 * xi + 128.0) * h1mxi;               
          d_eta_xi_d_xi = -65.0 / 128.0 * pow(xi - 1.0, 8.0) * (7.0 + xi * (56.0 + 3.0 * xi * (62.0 + xi * (104.0 + 77.0 * xi )))) * h1mxi;
          break;
        }
        default:
          io::messages.add("charge shape not implemented", "Lattice Sum",
                  io::message::critical);
      }
    }
    /**
     * calculates the fourier coefficient @f$ \hat{\gamma}(\kappa) @f$ of the 
     * charge shaping function
     * @param[in] shape charge shape switch
     * @param[in] kappa the value for which the coefficient is calculated
     * @param[out] gamma_hat the fourrer coefficient @f$ \hat{\gamma}(\kappa) @f$ 
     */
    static inline void charge_shape_fourier(const int &shape,
            const double &kappa,
            double &gamma_hat) {
      switch(shape) {
        case -1 :
          gamma_hat = exp(kappa * kappa * -0.25); 
          break;
        case 0 :
          gamma_hat = 3.0 * (sin(kappa) - kappa * cos(kappa))/pow(kappa,3.0);
          break;
        case 1 :
          gamma_hat = 12.0 * (2.0 - 2.0 * cos(kappa) - kappa * sin(kappa)) / pow(kappa,4.0);
          break;
        case 2 :
          gamma_hat = 60.0 * (2.0 * kappa + kappa * cos(kappa) - 3.0 * sin(kappa)) / pow(kappa,5.0);
          break;
        case 3 : 
          gamma_hat = 90.0*(8.0  + (kappa*kappa - 8.0)*cos(kappa) - 5.0 * kappa * sin(kappa))/pow(kappa,6.0);
          break;        
        case 4 :
          gamma_hat = 630.0 * (8.0 * kappa + 7.0 * kappa * cos(kappa) + (kappa * kappa - 15.0) * sin(kappa)) / pow(kappa,7.0);
          break;
        case 5: {
          const double k2 = kappa * kappa;
          gamma_hat = 5040.0 * (4.0 * (k2 - 6.0) - (k2 - 24.0) * cos(kappa) + 9.0 * kappa * sin(kappa)) / pow(kappa, 8.0);
          break;
        }
        case 6 : {
          const double k2 = kappa * kappa;
          gamma_hat = 7560.0 * ( 48.0 * kappa - (k2 - 57.0) * kappa * cos(kappa) + 3.0 * (4.0 * k2 - 35.0) * sin(kappa)) / pow(kappa, 9.0);
          break;
        }
        case 7 : {
          const double k2 = kappa * kappa;
          gamma_hat = 75600.0 * (24.0 * (k2 - 8.0) - 3.0 * (5.0 * k2 - 64.0) * cos(kappa) - (k2 - 87.0) * kappa * sin(kappa)) / pow(kappa, 10.0);
          break;
        }
        case 8 : {
          const double k2 = kappa * kappa;
          gamma_hat = 831600.0 * (8.0 * kappa * (k2 - 24.0) + (k2 - 123.0)* kappa * cos(kappa) - (18.0*k2 - 315.0) * sin(kappa)) / pow(kappa, 11.0);
          break;
        }
        case 9 : {
          const double k2 = kappa * kappa;
          gamma_hat = 1247400.0 * (192.0 * (k2 - 10.0) + (k2*k2 - 207.0*k2 + 1920) * cos(kappa) - (22.0 * k2 - 975.0) * kappa * sin(kappa)) / pow(kappa, 12.0);
          break;
        }
        case 10 : {
          const double k2 = kappa * kappa;
          gamma_hat = 16216200.0 * (64.0 * kappa * (k2 - 30.0) + (26.0 * k2 - 1545.0) * kappa * cos(kappa) + (k2*k2 - 285.0 * k2 + 3465.0) * sin(kappa)) / pow(kappa,13.0);
          break;
        }
        default:
          io::messages.add("charge shape not implemented", "Lattice Sum",
                  io::message::critical);
      }
    }    
    
    /**
     * one dimensinal charge assignment function @f$ w_p(\xi) @f$ 
     * @param[in] p order of the assignment function
     * @param[in] xi coordinate in grid units
     * @return the one-dimensional assignment function of order p
     */
    static inline double charge_assignment_1d(const int &p, const double &xi) {
      double absf = fabs(xi);
      switch (p) {
        case 1:
          return absf < 0.5 ? 1 : 0;
        case 2:
        {
          return absf < 1.0 ? (1.0 - absf) : 0;
        }
        case 3:
        {
          double result = 0.0;
          if (absf < 0.5)
            result = 0.75 - xi * xi;
          else if (absf >= 0.5 && absf < 1.5)
            result = 0.125 * pow(2.0 * absf - 3.0, 2.0);

          return result;
        }
        case 4:
        {
          double result = 0.0;
          if (absf < 1.0)
            result = 1.0 / 6.0 * (3.0 * pow(absf, 3.0) - 6.0 * pow(absf, 2.0) + 4.0);
          else if (absf >= 1.0 && absf < 2.0)
            result = -1.0 / 6.0 * pow(absf - 2.0, 3.0);

          return result;
        }
        case 5:
        {
          double result = 0.0;
          if (absf < 0.5)
            result = 1.0 / 192.0 * (48.0 * pow(xi, 4.0) - 120.0 * xi * xi + 115);
          else if (absf >= 0.5 && absf < 1.5)
            result = -1.0 / 96.0 * (16.0 * pow(xi, 4.0) - 80.0 * pow(absf, 3.0) +
                  120.0 * xi * xi - 20.0 * absf - 55.0);
          else if (absf >= 1.5 && absf < 2.5)
            result = 1.0 / 384.0 * pow(2 * absf - 5.0, 4.0);

          return result;
        }
        default:
          io::messages.add("assignment function not implemented", "Lattice Sum",
                  io::message::critical);
          return 0;
      }
    }
    
    /**
     * 3 dimensional charge assignment function @f$P_p(\mathbf{r}) = w_p(r_x) w_p(r_y) w_p(r_z)@f$
     * see MD.05.32 eq. 78
     * @param[in] p order of the assignment function
     * @param[in] a position vector in grid units
     * @return  the charge assignment function of the given point of order p 
     * 
     */
    static inline double charge_assignment(const int &p, const math::Vec & po) {
      return charge_assignment_1d(p, po(0)) * charge_assignment_1d(p, po(1)) * 
              charge_assignment_1d(p, po(2));
      
    }
    
    static inline double fourier_charge_assignment_1d(const int &p, const double &xi) {
      // MD.05.32 eq. 80
      if (xi == 0)
        return 1.0;
      else
        return pow(2.0 / xi * sin(xi / 2.0), p);
    }
    
    /**
     * returns the c parameters for a q as described in
     * MD99.32 p.65
     *
     * @param[in] q the order of the finite difference operator
     * @result a vector containing q finite difference constants.
     */
    static inline std::vector<double> finite_difference_constants(const int q) {
      std::vector<double> result(q, 0.0);
      switch (q) {
        case 0:
          break;
        case 1:
          result[0] = 1.0;
          break;
        case 2:
          result[0] = 4.0 / 3.0;
          result[1] = -1.0 / 3.0;
          break;
        case 3:
          result[0] = 3.0 / 2.0;
          result[1] = -3.0 / 5.0;
          result[2] = 0.1;
          break;
        case 4:
          result[0] = 8.0 / 5.0;
          result[1] = -4.0 / 5.0;
          result[2] = 8.0 / 35.0;
          result[3] = -1.0 / 35.0;
          break;
        case 5:
          result[0] = 5.0 / 3.0;
          result[1] = -20.0 / 21.0;
          result[2] = 5.0 / 14.0;
          result[3] = -5.0 / 63.0;
          result[4] = 1 / 126.0;
          break;
        default:
          io::messages.add("assignment function (difference constants) not implemented", "Lattice Sum",
                  io::message::critical);
      }
      return result;
    }
    
    /**
     * assign charges to grid
     *
     * The charges are taken from the topology and assigned to the charge 
     * densiy grid in the configutation using a charge assignment function of 
     * a given order.
     *
     * @param[in] r a VArray containing the positive, in-box atom positions
     */
    static void calculate_charge_density(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const math::VArray & r);


    /**
     * calculate the fourier transform of the electrostatic potential 
     * and the k-space energy
     *
     * From the charge density grid and the optimal influence function grid
     * in the configuration the electrostatic k-space energy is calculated. In 
     * the same loop the fourier transform of the potential is calculated 
     * and stored in the configuration.
     *
     * @param[in] storage the storage is used to store the energy
     */
    static void calculate_potential_and_energy(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            Storage & storage);
    
    /**
     * calculate the electric field
     *
     * From the fourier transform of the electrostatic potential the electric field
     * is calculated by ik or finite differentiation. The electrostatic field is
     * stored int the configuration
     */
    static void calculate_electric_field(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim);
    
    /**
     * calculate the force
     *
     * This function takes the electrostatic field from the configuration
     * and calculates the force using a similar assignment scheme as used in the
     * charge assignment procedure. 
     *
     * @param[in] storage the storage is used to store the force
     * @param[in] a VArray containing the positive, in-box positions of the charges
     */
    static void calculate_force(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            Storage & storage,
            const math::VArray & r);
  };
}

#endif

