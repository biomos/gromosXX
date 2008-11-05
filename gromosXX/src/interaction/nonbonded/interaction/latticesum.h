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
     * @todo optimize for speed 
     */
    static inline void charge_shape_switch(const int &shape,
            const double &xi,
            double &eta_xi,
            double &d_eta_xi_d_xi) {
      switch(shape) {
        case -1 :
          eta_xi = erfc(xi);
          d_eta_xi_d_xi = exp(-xi * xi) * -2.0 / sqrt(math::Pi);
          break;
        case 0 : {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          eta_xi = 0.5 * onemxi2*(xi + 2.0);
          d_eta_xi_d_xi = 3.0 / 2.0 * (xi * xi - 1.0);
          break;
        }
        case 1 : {
          const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi3 = onemxi2 * onemxi;
          eta_xi = onemxi3 * (xi + 1.0) ;
          d_eta_xi_d_xi = -2.0 * onemxi2 * (1.0 + 2.0 * xi);
          break;
        } 
        case 2 : {
          const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi3 = onemxi2 * onemxi;
          const double onemxi4 = onemxi2 * onemxi2;
          eta_xi = 1.0 / 2.0 * onemxi4 * (3.0 * xi + 2.0);
          d_eta_xi_d_xi = -5.0 / 2.0 * onemxi3 * (1.0 + 3.0 * xi);
          break;
        }
        case 3 : {
          const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi3 = onemxi2 * onemxi;
          const double onemxi4 = onemxi2 * onemxi2;
          eta_xi = 1.0 / 4.0 * onemxi4 * ((4.0 * xi + 7.0) * xi + 4.0);
          d_eta_xi_d_xi = -3.0 / 4.0 * onemxi3 * ((8.0 * xi + 9.0) * xi + 3.0);
          break;
        }
        case 4 : {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi4 = onemxi2 * onemxi2;
          const double onemxi5 = onemxi4 * onemxi;
          eta_xi = 1.0 / 8.0 * onemxi5 * ((15.0 * xi + 19.0) * xi + 8.0);
          d_eta_xi_d_xi = - 21.0 / 8.0 * onemxi4 * (1.0 + xi * ( 4.0 + 5.0 * xi));
          break;
        }
        case 5 : {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi4 = onemxi2 * onemxi2;
          const double onemxi5 = onemxi4 * onemxi;
          const double onemxi6 = onemxi4 * onemxi2;
          eta_xi = onemxi6 * ((3.0 * xi + 3.0) * xi + 1.0);
          d_eta_xi_d_xi = -3.0 * onemxi5 * (1.0 + xi * (5.0 + 8.0 * xi));
          break;
        }
        case 6 : {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi5 = onemxi2 * onemxi2 * onemxi;
          const double onemxi6 = onemxi5 * onemxi;
          eta_xi = 1.0 / 16.0 * onemxi6 * (((35.0 * xi + 66.0) * xi + 51.0) * xi + 16.0);
          d_eta_xi_d_xi = -9.0 / 16.0 * onemxi5 * (5.0 + xi * ( 25.0 + xi * (47.0 + 35.0 * xi )));
          break;
        }
        case 7 : {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi6 = onemxi2 * onemxi2 * onemxi2;
          const double onemxi7 = onemxi6 * onemxi;
          eta_xi = 1.0 / 8.0 * onemxi7 * (((32.0 * xi + 49.0) * xi + 31.0) * xi + 8.0);
          d_eta_xi_d_xi = -5.0 / 8.0 * onemxi6 * (5.0 + xi * (30.0 + xi * (69.0 + 64.0 * xi)));
          break;
        }
        case 8 : {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi3 = onemxi * onemxi * onemxi;
          const double onemxi7 = onemxi3 * onemxi3 * onemxi;
          const double onemxi8 = onemxi7 * onemxi;
          eta_xi = 1.0 / 16.0 * onemxi8 * (((105.0 * xi + 136.0) * xi + 73.0) * xi + 16.0);
          d_eta_xi_d_xi = -55.0 / 16.0 * onemxi7 * (1.0 + 3.0 * xi) * (1.0 + xi  * (4.0 + 7.0 * xi ));
          break;
        }
        case 9 : {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi3 = onemxi * onemxi * onemxi;
          const double onemxi7 = onemxi3 * onemxi3 * onemxi;
          const double onemxi8 = onemxi7 * onemxi;
          eta_xi = 1.0 / 32.0 * onemxi8 * ((((160.0 * xi + 335.0) * xi + 312.0) * xi + 151.0) * xi + 32.0);
          d_eta_xi_d_xi = -15.0 / 32.0 * onemxi7 * (7.0 + xi * (49.0 + xi * (141.0 + xi * (203.0 + 128.0*xi))));
          break;
        }
        case 10: {
         const double onemxi = 1.0 - xi;
          if (onemxi < 0.0) {
            eta_xi = 0.0;
            d_eta_xi_d_xi = 0.0;
            break;         
          }
          const double onemxi2 = onemxi * onemxi;
          const double onemxi4 = onemxi2 * onemxi2;
          const double onemxi8 = onemxi4 * onemxi4;
          const double onemxi9 = onemxi8 * onemxi;
          eta_xi =  1.0 / 128.0 * onemxi9 * ((((1155.0 * xi + 2075.0) * xi + 1665.0) * xi + 697.0) * xi + 128.0);               
          d_eta_xi_d_xi = -65.0 / 128.0 * onemxi8 * (7.0 + xi * (56.0 + 3.0 * xi * (62.0 + xi * (104.0 + 77.0 * xi ))));
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
        case -1 : {
          gamma_hat = exp(kappa * kappa * -0.25); 
          break;
        }
        case 0 : {
          const double k3 = kappa * kappa * kappa;
          gamma_hat = 3.0 * (sin(kappa) - kappa * cos(kappa))/k3;
          break;
        }
        case 1 : {
          const double k2 = kappa * kappa;
          const double k4 = k2*k2;
          gamma_hat = 12.0 * (2.0 - 2.0 * cos(kappa) - kappa * sin(kappa)) / k4;
          break;
        }
        case 2 : {
          const double k2 = kappa * kappa;
          const double k5 = k2 * k2 * kappa;
          gamma_hat = 60.0 * (2.0 * k2 * cos(kappa) - 3.0 * sin(kappa)) / k5;
          break;
        }
        case 3 : {
          const double k2 = kappa * kappa;
          const double k6 = k2 * k2 * k2;
          gamma_hat = 90.0*(8.0  + (k2 - 8.0)*cos(kappa) - 5.0 * kappa * sin(kappa))/k6;
          break;  
        }
        case 4 : {
          const double k2 = kappa;
          const double k7 = k2 * k2 * k2 * kappa;
          gamma_hat = 630.0 * (8.0 * kappa + 7.0 * kappa * cos(kappa) + (k2 - 15.0) * sin(kappa)) / k7;
          break;
        }
        case 5: {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k8 = k4 * k4;
          gamma_hat = 5040.0 * (4.0 * (k2 - 6.0) - (k2 - 24.0) * cos(kappa) + 9.0 * kappa * sin(kappa)) / k8;
          break;
        }
        case 6 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k9 = k4 * k4 * kappa;
          gamma_hat = 7560.0 * ( 48.0 * kappa - (k2 - 57.0) * kappa * cos(kappa) + 3.0 * (4.0 * k2 - 35.0) * sin(kappa)) / k9;
          break;
        }
        case 7 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k10 = k4 * k4 * k2;
          gamma_hat = 75600.0 * (24.0 * (k2 - 8.0) - 3.0 * (5.0 * k2 - 64.0) * cos(kappa) - (k2 - 87.0) * kappa * sin(kappa)) / k10;
          break;
        }
        case 8 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k11 = k4 * k4 * k2 * kappa;
          gamma_hat = 831600.0 * (8.0 * kappa * (k2 - 24.0) + (k2 - 123.0)* kappa * cos(kappa) - (18.0*k2 - 315.0) * sin(kappa)) / k11;
          break;
        }
        case 9 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k12 = k4 * k4 * k4;
          gamma_hat = 1247400.0 * (192.0 * (k2 - 10.0) + (k4 - 207.0*k2 + 1920) * cos(kappa) - (22.0 * k2 - 975.0) * kappa * sin(kappa)) / k12;
          break;
        }
        case 10 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k13 = k4 * k4 * k4 * kappa;
          gamma_hat = 16216200.0 * (64.0 * kappa * (k2 - 30.0) + (26.0 * k2 - 1545.0) * kappa * cos(kappa) + (k4 - 285.0 * k2 + 3465.0) * sin(kappa)) / k13;
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
          if (absf < 0.5) {
            return 0.75 - xi * xi;
          } else if (absf >= 0.5 && absf < 1.5) {
            const double term = 2.0 * absf - 3.0;
            return 0.125 * term * term;
          } else {
            return 0.0;
          }
        }
        case 4:
        {
          if (absf < 1.0) {
            const double absf2 = absf * absf;
            return 1.0 / 6.0 * (3.0 * absf2 * absf - 6.0 * absf2 + 4.0);
          } else if (absf >= 1.0 && absf < 2.0) {
            const double term = absf - 2.0;
            return -1.0 / 6.0 * term * term * term;
          } else {
            return 0.0;
          }
        }
        case 5:
        {
          if (absf < 0.5) {
            const double xi2 = xi * xi;
            const double xi4 = xi2 * xi2;
            return 1.0 / 192.0 * (48.0 * xi4 - 120.0 * xi2 + 115);
          } else if (absf >= 0.5 && absf < 1.5) {
            const double xi2 = xi * xi;
            const double xi4 = xi2 * xi2;
            const double absf3 = absf * absf * absf;
            return -1.0 / 96.0 * (16.0 * xi4 - 80.0 * absf3 +
                  120.0 * xi2 - 20.0 * absf - 55.0);
          } else if (absf >= 1.5 && absf < 2.5) {
            const double tmp = 2 * absf - 5.0;
            const double tmp2 = tmp * tmp;
            return 1.0 / 384.0 * tmp2 * tmp2;
          } else {
            return 0.0;
          }
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
     * the @f$f_p@f$ function for the P3M self term evaluation
     */
    static inline double p3m_selfterm_fp(int p, int n, double s) {
      switch (p) {
        case 1:
        {
          return 1.0;
        } // p = 1
        case 2:
        {
          const double s2 = s*s;
          switch (abs(n)) {
            case 0: return 1.0 / 2.0 * (1.0 + 4.0 * s2);
            case 1: return 1.0 / 4.0 * (1.0 - 4.0 * s2);
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp", io::message::critical);
              return 0.0;
          }
        } // p = 2
        case 3:
        {
          const double s2 = s*s;
          switch (abs(n)) {
            case 0: return 1.0 / 32.0 * (19.0 + s2 * (-24.0 + s2 * 48.0));
            case 1: return 1.0 / 16.0 * (3.0 + s2 * (8.0 - s2 * 16.0));
            case 2:
            {
              const double term = 1.0 - 4.0 * s2;
              return 1.0 / 64.0 * term *term;
            }
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp", io::message::critical);
              return 0.0;
          }
        }// p = 3
        case 4:
        {
          const double s2 = s*s;
          switch (abs(n)) {
            case 0: return 1.0 / 576.0 * (265.0 + s2 * (204.0 + s2 *
                      (-528.0 + s2 * 320.0)));
            case 1: return 1.0 / 2304.0 * (575.0 + s2 * (-564.0 + s2 *
                      (1488.0 + s2 * -960.0)));
            case 2: return 1.0 / 1152.0 * (23.0 + s2 * (84.0 + s2 *
                      (-240.0 + s2 * 192.0)));
            case 3:
            {
              const double term = 1.0 - 4.0 * s2;
              return 1.0 / 2304.0 * term * term *term;
            }
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp", io::message::critical);
              return 0.0;
          }
        } // p = 4
        case 5:
        {
          const double s2 = s*s;
          switch (abs(n)) {
            case 0: return 1.0 / 73728.0 * (32227.0 + s2 * (-9520.0 + s2 *
                      (28960.0 + s2 * (-29440.0 + s2 * 8960.0))));
            case 1: return 1.0 / 18432.0 * (4389.0 + s2 * (1792.0 + s2 *
                      (-5472.0 + s2 * (5632.0 + s2 * -1792.0))));
            case 2: return 1.0 / 36864.0 * (1559.0 + s2 * (-1456.0 + s2 *
                      (4512.0 + s2 * (-4864.0 + s2 * 1792.0))));
            case 3: return 1.0 / 18432.0 * (19.0 + s2 * (128.0 + s2 *
                      (-416.0 + s2 * (512.0 + s2 * -256.0))));
            case 4:
            {
              const double term = 1.0 - 4.0 * s2;
              const double term2 = term * term;
              return 1.0 / 147456.0 * term2 * term2;
            }
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp", io::message::critical);
              return 0.0;
          }
        } // p = 5
          io::messages.add("fp: p invalid.", "p3m_selfterm_fp", io::message::critical);
          return 0.0;
      }
      return 0.0;
    }
    
     /**
     * the averaged @f$f_p@f$ function for the P3M self term evaluation
     */
    static inline double p3m_selfterm_fp_avg(int p, int n) {
      switch (p) {
        case 1:
        {
          return 1.0;
        } // p = 1
        case 2:
        {
          switch (abs(n)) {
            case 0: return 2.0 / 3.0;
            case 1: return 1.0 / 6.0;
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp_avg", io::message::critical);
              return 0.0;
          }
        } // p = 2
        case 3:
        {
          switch (abs(n)) {
            case 0: return 11.0 / 20.0;
            case 1: return 13.0 / 60.0;
            case 2: return 1.0 / 120.0;
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp_avg", io::message::critical);
              return 0.0;
          }
        } // p = 3
        case 4:
        {
          switch (abs(n)) {
            case 0: return 151.0 / 315.0;
            case 1: return 397.0 / 1680.0;
            case 2: return 1.0 / 42.0;
            case 3: return 1.0 / 5040.0;
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp_avg", io::message::critical);
              return 0.0;
          }
        } // p = 4
        case 5:
        {
          switch (abs(n)) {
            case 0: return 15619.0 / 36288.0;
            case 1: return 44117.0 / 181440.0;
            case 2: return 913.0 / 22680.0;
            case 3: return 251.0 / 181440.0;
            case 4: return 1.0 / 362880.0;
            default:
              io::messages.add("fp: |n| invalid.", "p3m_selfterm_fp_avg", io::message::critical);
              return 0.0;
          }
        } // p = 5
      }
      return 0.0;
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
    template<class MeshType>
    static void calculate_charge_density(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const math::VArray & r);
    
    /**
     * assign square charges to grid
     *
     * The squared charges are taken from the topology and assigned to the  
     * @f$R_g@f$ grid in the configutation using a charge assignment function of 
     * a given order.
     *
     * @param[in] r a VArray containing the positive, in-box atom positions
     */
    template<class MeshType>
    static void calculate_squared_charge_grid(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const math::VArray & r);
    
     /**
     * assign square charges to grid
     *
     * The squared charges are taken from the topology and assigned to the  
     * @f$R_g@f$ grid in the configutation using an averaged charge assignment
     * function of a given order.
     */
    template<class MeshType>
    static void calculate_averaged_squared_charge_grid(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim);


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
    template<class MeshType>
    static void calculate_potential_and_energy(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            Storage & storage);
    
    /**
     * calculates the p3m ~A2 term
     *
     * @param[in] storage the storage is used to store the energy
     */
    template<class MeshType>
    static void calculate_p3m_selfterm(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim);
    
    /**
     * calculate the electric field
     *
     * From the fourier transform of the electrostatic potential the electric field
     * is calculated by ik or finite differentiation. The electrostatic field is
     * stored int the configuration
     */
    template<class MeshType>
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
    template<class MeshType>
    static void calculate_force(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            Storage & storage,
            const math::VArray & r);
    
    /**
     * decompose domains
     *
     * For every atom calculate the x component of the grid assignment. If it
     * lies in the current domain, add the atom the domain.
     */
    template<class MeshType>
    static void decompose_into_domains(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const math::VArray & r,
            const unsigned int size);
  };
 
  
}

#endif

