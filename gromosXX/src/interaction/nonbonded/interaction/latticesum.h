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
      Lattice_Sum();
  public:     
    /**
     * the virtual function to ensure the class is not instantiated.
     */
    //virtual void dummy() = 0;
    
    /**
     * calculates charge shaping switch function @f$ \eta(\xi) @f$ and its derivative
     * @f$ \eta'(\xi) @f$ 
     *
     * depending on the shape the following functions are evalulated
     * <table border="0">
     * <tr>
     *   <td><em>shape</em> @f$ N_{\gamma} @f$</td>
     *   <td><em>switch function</em> @f$ \eta(\xi) @f$</td>
     *   <td><em>derivative</em> @f$ -\eta'(\xi) @f$</td></tr>
     * <tr><td>-1</td><td> @f$\mathrm{erfc}(\xi) @f$ </td>
     *   <td>@f$ 2\pi^{-1/2}e^{-\xi^2} @f$</td></tr>
     * <tr><td>0</td><td> @f$ (1/2)\,(1-\xi)^2(\xi+2)\,H(1-\xi) @f$ </td>
     *   <td>@f$ (3/2)\,(1-\xi)(\xi+1)\,H(1-\xi) @f$</td></tr>
     * <tr><td>1</td><td> @f$ (1-\xi)^3(\xi+1)\,H(1-\xi)@f$ </td>
     *   <td>@f$ 2(1-\xi)^2(2\xi+1)\,H(1-\xi) @f$</td></tr>
     * <tr><td>2</td><td> @f$ (1/2)\,(1-\xi)^4(3\xi+2)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (5/2)\,(1-\xi)^3(3\xi+1)\,H(1-\xi) @f$</td></tr>
     * <tr><td>3</td><td> @f$ (1/4)\,(1-\xi)^4(4\xi^2+7\xi+4)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (3/4)\,(1-\xi)^3(8\xi^2+9\xi+3)\,H(1-\xi) @f$</td></tr>
     * <tr><td>4</td><td> @f$(1/8)(1-\xi)^5(15\xi^2+19\xi+8)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (21/8)(1-\xi)^4(5\xi^2+4\xi+1)\,H(1-\xi) @f$</td></tr>
     * <tr><td>5</td><td> @f$ (1-\xi)^6(3\xi^2+3\xi+1)\,H(1-\xi)@f$ </td>
     *   <td>@f$ 3(1-\xi)^5(8\xi^2+5\xi+1)\,H(1-\xi) @f$</td></tr>
     * <tr><td>6</td><td> @f$(1/16)\,(1-\xi)^6(35\xi^3+66\xi^2+51\xi+16)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (9/16)\,(1-\xi)^5(35\xi^3+47\xi^2+25\xi+5)\,H(1-\xi) @f$</td></tr>
     * <tr><td>7</td><td> @f$(1/8)(1-\xi)^7(32\xi^3+49\xi^2+31\xi+8)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (5/8)(1-\xi)^6(64\xi^3+69\xi^2+30\xi+5)\,H(1-\xi) @f$</td></tr>
     * <tr><td>8</td><td> @f$(1/16)\,(1-\xi)^8(105\xi^3+136\xi^2+73\xi+16)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (55/16)\,(1-\xi)^7(21\xi^3+19\xi^2+7\xi+1)\,H(1-\xi) @f$</td></tr>
     * <tr><td>9</td><td> @f$(1/32)\,(1-\xi)^8(160\xi^4+335\xi^3+312\xi^2+151\xi+32)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (15/32)\,(1-\xi)^7(128\xi^4+203\xi^3+141\xi^2+49\xi+7)\,H(1-\xi) @f$</td></tr>
     * <tr><td>10</td><td> @f$(1/128)(1-\xi)^9(1155\xi^4+2075\xi^3+1665\xi^2+697\xi+128)\,H(1-\xi)@f$ </td>
     *   <td>@f$ (65/128)(1-\xi)^8(231\xi^4+312\xi^3+186\xi^2+56\xi+7)\,H(1-\xi) @f$</td></tr>
     * </table>
     * where @f$H@f$ is the heaviside step function.
     *
     * @param[in] shape the charge shape switch
     * @param[out] eta_xi the shaping switch function @f$ \eta(\xi) @f$
     * @param[out] eta_xi_d_xi the derivative @f$ \eta'(\xi) @f$
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
     * charge shaping function. If a pointer for gamma_hat_prime is given,
     * the derivative @f$ \hat{\gamma}'(\kappa) @f$ is calculated as well.
     *
     * <table border="0">
     * <tr>
     *   <td><em>shape</em> @f$ N_{\gamma} @f$</td>
     *   <td><em>fourier coefficient</em> @f$ \kappa^{N_{\gamma}+3} \hat{\gamma}(\kappa) @f$ </td>
     *   <td><em>derivative</em> @f$ \kappa^{N_{\gamma}+4} \hat{\gamma}'(\kappa) @f$ </td>
     * </tr>
     * <tr><td>-1</td>
     *   <td>@f$ \kappa^2 e^{-\kappa^2/4} @f$</td>
     *   <td>@f$ -(1/2)\kappa^4  e^{-\kappa^2/4} @f$</td></tr>
     * <tr><td>0</td>
     *   <td>@f$ 3[-\kappa C+S] @f$</td>
     *   <td>@f$ 3[3\kappa C+(\kappa^2-3)S] @f$</td></tr>
     * <tr><td>1</td>
     *   <td>@f$ 12[2-2C-\kappa S] @f$</td>
     *   <td>@f$ 12[-8-(\kappa^2-8)C+5\kappa S] @f$</td></tr>
     * <tr><td>2</td>
     *   <td>@f$ 60[2\kappa +\kappa C-3S] @f$</td>
     *   <td>@f$ 60[-8\kappa - 7\kappa C - (\kappa^2-15)S] @f$</td></tr>
     * <tr><td>3</td>
     *   <td>@f$ 90[8+(\kappa ^2-8)C-5\kappa S] @f$</td>
     *   <td>@f$ 90[-48 - (9\kappa ^2-48)C - (\kappa^2-33)\kappa S] @f$</td></tr>
     * <tr><td>4</td>
     *   <td>@f$ 630[8\kappa +7\kappa C+(\kappa ^2-15)S] @f$</td>
     *   <td>@f$ 630[-48\kappa + (\kappa^2-57) \kappa C - 3(4\kappa ^2-35)S] @f$</td></tr>
     * <tr><td>5</td>
     *   <td>@f$ 5040[4(\kappa ^2-6)-(\kappa ^2-24)C+9\kappa S] @f$</td>
     *   <td>@f$ 5040[-24(\kappa ^2-8)+3(5\kappa ^2-64)C+(\kappa^2-87)\kappa S] @f$</td></tr>
     * <tr><td>6</td>
     *   <td>@f$ 7560[48\kappa -(\kappa ^2-57)\kappa C+3(4\kappa ^2-35)S] @f$</td>
     *   <td>@f$ 7560[-384\kappa +3(6\kappa^2-187)\kappa C+(\kappa^4-141\kappa^2+945)S] @f$</td></tr>
     * <tr><td>7</td>
     *   <td>@f$ 75600[24(\kappa ^2-8)-3(5\kappa ^2-64)C-(\kappa ^2-87)\kappa S] @f$</td>
     *   <td>@f$ 75600[-192(\kappa ^2-10)-(\kappa^4-207\kappa^2+1920)C+(22\kappa^2-975)\kappa S] @f$</td></tr>
     * <tr><td>8</td>
     *   <td>@f$ 831600[8\kappa (\kappa ^2-24)+(\kappa ^2-123)\kappa C-(18\kappa ^2-315)S] @f$</td>
     *   <td>@f$ 831600[-64\kappa (\kappa^2-30)-(26\kappa ^2-1545)\kappa C-(\kappa ^4-285\kappa^2+3465)S] @f$</td></tr>
     * <tr><td>9</td>
     *   <td>@f$ 1247400[192(\kappa ^2-10)+(\kappa ^4-207\kappa ^2+1920)C-(22\kappa ^2-975)\kappa S] @f$</td>
     *   <td>@f$ 1247400[-1920(\kappa ^2-12)-15(2\kappa ^4-203\kappa ^2+1536)C-(\kappa^4-405\kappa ^2+12645)\kappa S] @f$</td></tr>
     * <tr><td>10</td>
     *   <td>@f$ 16216200[64\kappa (\kappa ^2-30)+(26\kappa ^2-1545)\kappa C+(\kappa ^4-285\kappa^2+3465)S] @f$</td>
     *   <td>@f$ 16216200[-640\kappa (\kappa ^2-36)+(\kappa^4-545\kappa ^2+22005)\kappa C-5(7\kappa ^4-936\kappa^2+9009)S] @f$</td></tr>
     * </table>
     * with @f$S = \sin(\kappa)@f$ and @f$C = \cos(\kappa)@f$
     *
     * @param[in] shape charge shape switch
     * @param[in] kappa the value for which the coefficient is calculated
     * @param[out] gamma_hat the fourrer coefficient @f$ \hat{\gamma}(\kappa) @f$ 
     * @param[out] gamma_hat_prime the derivative of the fourrier coefficient
     *             (pointer as it can be omitted by giving NULL)
     */
    static inline void charge_shape_fourier(const int &shape,
            const double &kappa,
            double &gamma_hat, double * gamma_hat_prime = NULL) {
      switch(shape) {
        case -1 : {
          const double k2_quater = kappa * kappa * -0.25;
          gamma_hat = exp(k2_quater); 
          if (gamma_hat_prime != NULL) {
            *gamma_hat_prime = -0.5 * kappa * exp(k2_quater);
          }
          break;
        }
        case 0 : {
          const double k2 = kappa * kappa;
          const double k3 = k2 * kappa;
          const double S = sin(kappa);
          const double C = cos(kappa);
          const double kappa_C = kappa * C;
          gamma_hat = 3.0 * (S - kappa_C)/k3;
          if (gamma_hat_prime != NULL) {
            const double k4 = k2 * k2;
            *gamma_hat_prime = 3.0 * (3.0 * kappa_C + (k2 - 3.0)*S) / k4;
          }
          break;
        }
        case 1 : {
          const double k2 = kappa * kappa;
          const double k4 = k2*k2;
          const double C = cos(kappa);
          const double kappa_S = kappa*sin(kappa);
          gamma_hat = 12.0 * (2.0 - 2.0 * C - kappa_S) / k4;
          if (gamma_hat_prime != NULL) {
            const double k5 = k4 * kappa;
            *gamma_hat_prime = 12.0 * (-8.0 - (k2 - 8.0)*C + 5.0 * kappa_S) / k5;
          }
          break;
        }
        case 2 : {
          const double k2 = kappa * kappa;
          const double k5 = k2 * k2 * kappa;
          const double C = cos(kappa);
          const double S = sin(kappa);
    //   wrong!!!   gamma_hat = 60.0 * (2.0 * k2 * C - 3.0 * S) / k5;
          gamma_hat = 60.0 * (2.0 * kappa + kappa * C - 3.0 * S) / k5;
          if (gamma_hat_prime != NULL) {
            const double k6 = k5 * kappa;
            *gamma_hat_prime = 60.0 * (-8.0 * kappa - 7.0 * kappa * C - (k2 - 15.0)*S) / k6;
          }
          break;
        }
        case 3 : {
          const double k2 = kappa * kappa;
          const double k6 = k2 * k2 * k2;
          const double C = cos(kappa);
          const double kappa_S = kappa * sin(kappa);
          gamma_hat = 90.0*(8.0  + (k2 - 8.0)*C - 5.0 * kappa_S)/k6;
          if (gamma_hat_prime != NULL) {
            const double k7 = k6 * kappa;
            *gamma_hat_prime = 90.0 * (-48.0 - (9.0*k2 - 48.0)*C - (k2 - 33.0)*kappa_S)/k7;
          }
          break;  
        }
        case 4 : {
         // wrong!!! const double k2 = kappa ;
          const double k2 = kappa * kappa ;
          const double k7 = k2 * k2 * k2 * kappa;
          const double kappa_C = kappa * cos(kappa);
          const double S = sin(kappa);
          gamma_hat = 630.0 * (8.0 * kappa + 7.0 * kappa_C + (k2 - 15.0) * S) / k7;
          if (gamma_hat_prime != NULL) {
            const double k8 = k7 * kappa;
           // wrong!!! *gamma_hat_prime = 630.0 * (-48.0*kappa + (k2 - 57.0)*kappa_C - 3.0 * (4.0*k2 - 35.0)*kappa*S)/k8;
            *gamma_hat_prime = 630.0 * (-48.0*kappa + (k2 - 57.0)*kappa_C - 3.0 * (4.0*k2 - 35.0)*S)/k8;
          }
          break;
        }
        case 5: {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k8 = k4 * k4;
          const double C = cos(kappa);
          const double kappa_S = kappa*sin(kappa);
          gamma_hat = 5040.0 * (4.0 * (k2 - 6.0) - (k2 - 24.0) * C + 9.0 * kappa_S) / k8;
          if (gamma_hat_prime != NULL) {
            const double k9 = k8 * kappa;
            *gamma_hat_prime = 5040.0 * (-24.0 * (k2 - 8.0) + 3.0 * (5.0 * k2 - 64.0) * C + (k2 - 87.0) * kappa_S) / k9;
          }
          break;
        }
        case 6 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k9 = k4 * k4 * kappa;
          const double kappa_C = kappa * cos(kappa);
          const double S = sin(kappa);
          gamma_hat = 7560.0 * ( 48.0 * kappa - (k2 - 57.0) * kappa_C + 3.0 * (4.0 * k2 - 35.0) * S) / k9;
          if (gamma_hat_prime != NULL) {
            const double k10 = k9 * kappa;
            *gamma_hat_prime = 7560.0 * (-384.0 * kappa + 3.0 * (6.0 * k2 - 187.0) * kappa_C + (k4 - 141.0 * k2 + 945.0) * S) / k10;
          }
          break;
        }
        case 7 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k10 = k4 * k4 * k2;
          const double C = cos(kappa);
          const double kappa_S = kappa * sin(kappa);
          gamma_hat = 75600.0 * (24.0 * (k2 - 8.0) - 3.0 * (5.0 * k2 - 64.0) * C - (k2 - 87.0) * kappa_S) / k10;
          if (gamma_hat_prime != NULL) {
            const double k11 = k10 * kappa;
            *gamma_hat_prime = 75600.0 * (-192.0 * (k2 - 10.0)-(k4 - 207.0 * k2 + 1920.0) * C + (22.0 * k2 - 975.0) * kappa_S) / k11;
          }
          break;
        }
        case 8 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k11 = k4 * k4 * k2 * kappa;
          const double kappa_C = kappa * cos(kappa);
          const double S = sin(kappa);
          gamma_hat = 831600.0 * (8.0 * kappa * (k2 - 24.0) + (k2 - 123.0)* kappa_C - (18.0*k2 - 315.0) * S) / k11;
          if (gamma_hat_prime != NULL) {
            const double k12 = k11 * kappa;
            *gamma_hat_prime = 831600.0 * (-64.0 * kappa * (k2 - 30.0)-(26.0 * k2 - 1545.0) * kappa_C - (k4 - 285.0 * k2 + 3465.0) * S) / k12;
          }
          break;
        }
        case 9 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k12 = k4 * k4 * k4;
          const double C = cos(kappa);
          const double kappa_S = kappa * sin(kappa);
          gamma_hat = 1247400.0 * (192.0 * (k2 - 10.0) + (k4 - 207.0*k2 + 1920) * C - (22.0 * k2 - 975.0) * kappa_S) / k12;
          if (gamma_hat_prime != NULL) {
            const double k13 = k12 * kappa;
            *gamma_hat_prime = 1247400.0 * (-1920.0 * (k2 - 12.0) - 15.0 * (2.0 * k4 - 203.0 * k2 + 1536.0) * C - (k4 - 405.0 * k2 + 12645.0) * kappa_S) / k13;
          }
          break;
        }
        case 10 : {
          const double k2 = kappa * kappa;
          const double k4 = k2 * k2;
          const double k13 = k4 * k4 * k4 * kappa;
          const double kappa_C = kappa * cos(kappa);
          const double S = sin(kappa);
          gamma_hat = 16216200.0 * (64.0 * kappa * (k2 - 30.0) + (26.0 * k2 - 1545.0) * kappa_C + (k4 - 285.0 * k2 + 3465.0) * S) / k13;
          if (gamma_hat_prime != NULL) {
            const double k14 = k13 * kappa;
            *gamma_hat_prime = 16216200.0 * (-640.0 * kappa * (k2 - 36.0)+(k4 - 545.0 * k2 + 22005.0) * kappa_C - 5.0 * (7.0 * k4 - 936.0 * k2 + 9009.0) * S) / k14;
          }
          break;
        }
        default:
          io::messages.add("charge shape not implemented", "Lattice Sum",
                  io::message::critical);
      }
    }
    
    /**
     * one dimensinal charge assignment function @f$ w_p(\xi) @f$
     * <table border="0">
     * <tr><td><em>order</em> @f$ p @f$</td>
     *   <td><em>assignment function</em> @f$ w_p(\xi) @f$</td>
     *   <td><em>range</em></td></tr>
     * <tr><td>1</td>
     *   <td>@f$ 1 @f$</td>
     *   <td>@f$ |\xi| < 1/2 @f$</td></tr>
     * <tr><td>2</td>
     *   <td>@f$ 1-|\xi| @f$</td>
     *   <td>@f$ |\xi| < 1 @f$</td></tr>
     * <tr><td>3</td>
     *   <td>@f$ -(3/4)(4\xi^2-1) @f$</td>
     *   <td>@f$ |\xi| < 1/2 @f$</td></tr>
     * <tr><td>&nbsp;</td>
     *   <td>@f$ (1/8)(2|\xi|-3)^2 @f$</td>
     *   <td>@f$ 1/2 < |\xi| < 3/2 @f$</td></tr>
     * <tr><td>4</td>
     *   <td>@f$ (1/6)(3|\xi|^3-6 \xi^2+4) @f$</td>
     *   <td>@f$ |\xi| < 1 @f$</td></tr>
     * <tr><td>&nbsp;</td>
     *   <td>@f$ -(1/6)(|\xi|-2)^3 @f$</td>
     *   <td>@f$ 1 < |\xi| < 2 @f$</td></tr>
     * <tr><td>5</td>
     *   <td>@f$ (1/192)(48\xi^4-120 \xi^2+115) @f$</td>
     *   <td>@f$ |\xi| < 1/2 @f$</td></tr>
     * <tr><td></td>
     *   <td>@f$ -(1/96)(16\xi^4-80|\xi|^3+120\xi^2-20|\xi|-55) @f$</td>
     *   <td>@f$ 1/2 < |\xi| < 3/2 @f$</td></tr>
     * <tr><td></td>
     *   <td>@f$ (1/384)(2|\xi|-5)^4 @f$</td>
     *   <td>@f$ 3/2 < |\xi| < 5/2 @f$</td></tr>
     * </table>
     * 
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
     * with @f$ w_p @f$ as  one-dimensional charge assignment function.
     *
     * @param[in] p order of the assignment function
     * @param[in] a position vector in grid units
     * @return  the charge assignment function of the given point of order p
     * @sa MD.05.32 eq. 78 @ref interaction::Lattice_Sum::charge_assignment_1d
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
     *
     * <table border="0">
     *   <tr><td><em>order</em> @f$ p @f$</td><td>@f$ n @f$ </td>
     *    <td>@f$ f_p(n; s) @f$</td></tr>
     *   <tr><td>1</td><td>0</td><td>@f$ 1 @f$</td></tr>
     *   <tr><td>2</td><td>0</td><td>@f$ 1/2(1 + 4s^2) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>1</td><td>@f$ 1/4(1 - 4s^2) @f$</td></tr>
     *   <tr><td>3</td><td>0</td><td>@f$ 1/32(19-24 s^2 + 48s^4) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>1</td><td>@f$ 1/16(3+8s^2-16s^4) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>2</td><td>@f$ 1/64(1-8s^2+16s^4)=1/64(1-4s^2)^2 @f$</td></tr>
     *   <tr><td>4</td><td>0</td><td>@f$ 1/576(265+204s^2-528s^4+320s^6) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>2</td><td>@f$ 1/2304(575-564s^2+1488s^4-960s^6) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>3</td><td>@f$ 1/1152(23+84s^2-240s^4+192s^6) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>4</td><td>@f$ 1/2304(1-12s^2+48s^4-64s^6)=1/2304(1-4s^2)^3 @f$</td></tr>
     *   <tr><td>5</td><td>0</td><td>@f$ 1/73728 (32227-9520s^2+28960s^4-29440s^6+8960s^8) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>1</td><td>@f$ 1/18432(4389+1792s^2-5472s^4+5632s^6-1792s^8) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>2</td><td>@f$ 1/36864(1559-1456s^2+4512s^4-4864s^6+1792s^8) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>3</td><td>@f$ 1/18432(19+128s^2-416s^4+512s^6-256s^8) @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>4</td><td>@f$ 1/147456(1-16s^2+96s^4-256s^6+256s^8)=1/147456(1-4s^2)^4 @f$</td></tr>
     * </table>
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
        default :
        {
          io::messages.add("fp: p invalid.", "p3m_selfterm_fp", io::message::critical);
        }
      }
      return 0.0;
    }

    /**
     * the averaged @f$f_p@f$ function for the P3M self term evaluation
     *
     * <table border="0">
     *   <tr><td><em>order</em> @f$ p @f$</td><td>@f$ n @f$ </td>
     *    <td>@f$ \overline{f}_p(n) @f$</td></tr>
     *   <tr><td>1</td><td>0</td><td>@f$ 1 @f$</td></tr>
     *   <tr><td>2</td><td>0</td><td>@f$ 2/3 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>1</td><td>@f$ 1/6 @f$</td></tr>
     *   <tr><td>3</td><td>0</td><td>@f$ 11/20 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>1</td><td>@f$ 13/60 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>2</td><td>@f$ 1/120 @f$</td></tr>
     *   <tr><td>4</td><td>0</td><td>@f$ 151/315 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>2</td><td>@f$ 397/1680 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>3</td><td>@f$ 1/42 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>4</td><td>@f$ 1/5040 @f$</td></tr>
     *   <tr><td>5</td><td>0</td><td>@f$ 15619/36288 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>1</td><td>@f$ 44117/181440 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>2</td><td>@f$ 913/22680 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>3</td><td>@f$ 251/181440 @f$</td></tr>
     *   <tr><td>&nbsp;</td><td>4</td><td>@f$ 1/362880 @f$</td></tr>
     * </table>
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
     * The charge density @f$ s_g @f$ grid is calculated as
     * @f[ s_g(\mathbf{r}_n) = \sum_{i=1}^{N_q} P(\mathbf{r}_n - \mathbf{r}_i) @f]
     *
     * where @f$ P(\mathbf{r}) @f$ is the @ref interaction::Lattice_Sum::charge_assignment 
     * "assignment function"
     *
     * @param[in] domain the indices of the atoms to use.
     * @param[in] r a VArray containing the positive, in-box atom positions
     */
    template<class MeshType>
    static void calculate_charge_density(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const std::vector<int> & domain,
            const math::VArray & r);
    
    /**
     * assign square charges to grid
     *
     * The squared charges are taken from the topology and assigned to the  
     * @f$ R_g @f$ grid in the configutation using a charge assignment function of 
     * a given order.
     *
     * @f[ R_g(\mathbf{r}_{\mathbf{n}}) = V_G^{-1} \sum_{i=1}^{N_q} q_{i}^2
              \sum_{\mathbf{m} \in Z^3} \prod_{\mu=x,y,z}
              f_p(L_{\mu}m_{\mu}+n_\mu;s_{i,\mu}) @f]
     * with @f$ \mathbf{s}_i = \underline{H}^{-1}\mathbf{r}_i @f$ being the 
     * distance of the charge @f$ i @f$ to the closest grid point or centre. This
     * assigns the squared charge to @f$ (2p-1)^3 @f$ grid points neighboring the
     * origin of the grid.
     *
     * For the function @f$ f_p(n; s) @f$ see 
     * @ref interaction::Lattice_Sum::p3m_selfterm_fp "here" .
     *
     * @param[in] domain the indices of the atoms to use
     * @param[in] r a VArray containing the positive, in-box atom positions
     * @sa interaction::Lattice_Sum::p3m_selfterm_fp 
     * @sa interaction::Lattice_Sum::calculate_averaged_squared_charge_grid
     */
    template<class MeshType>
    static void calculate_squared_charge_grid(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const std::vector<int> & domain,
            const math::VArray & r);
    
    /**
     * assign averaged square charges to grid
     *
     * The squared charges are taken from the topology and assigned to the  
     * @f$R_g@f$ grid in the configutation using an averaged charge assignment
     * function of a given order.
     *
     * See @ref interaction::Lattice_Sum::calculate_squared_charge_grid for 
     * details. The difference is, that @f$ \overline{f}_p(n) @f$ is used
     * as assignment function.
     *
     * @param[in] domain the indices of thea toms to use
     * @sa interaction::Lattice_Sum::calculate_squared_charge_grid
     * @sa interaction::Lattice_Sum::p3m_selfterm_fp_avg
     */
    template<class MeshType>
    static void calculate_averaged_squared_charge_grid(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const std::vector<int> & domain);


    /**
     * calculate the fourier transform of the electrostatic potential 
     * and the k-space energy
     *
     * From the charge density grid and the optimal influence function grid
     * in the configuration the electrostatic k-space energy is calculated. In 
     * the same loop the fourier transform of the potential is calculated 
     * and stored in the configuration.
     *
     * If required, the virial contribution is calculated.
     *
     * The energy is calculated as
     * @f[ E_{\gamma}  = (2\epsilon_0V)^{-1} \sum_{\mathbf{l} \in G, l\neq 0}
             \hat{G}_g^{\dagger}(\mathbf{k}_{\mathbf{l}})
             |\hat{s}_g(\mathbf{k}_{\mathbf{l}})|^2 @f]
     *
     * The fourier transform of the potential is calculated as
     * @f[ \hat{\Phi}_{\gamma,g}(\mathbf{k}_{\mathbf{l}}) = \epsilon_0^{-1}
             \hat{G}_g^{\dagger}(\mathbf{k}_{\mathbf{l}})
             \hat{s}_g(\mathbf{k}_{\mathbf{l}}) @f]
     *
     * The virial is calculated as
     * @f[ \underline{W}_\gamma = -(4\epsilon_0V)^{-1} \sum_{\mathbf{l} \in G, l\neq 0}
             \left[\hat{\underline{\Gamma}}_g^0(\mathbf{k}_{\mathbf{l}}) + 
             \mathrm{diag}(\hat{G}_g^{\dagger}(\mathbf{k}_{\mathbf{l}}))\right] 
             |\hat{s}_g(\mathbf{k}_{\mathbf{l}})|^2 @f]   
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
     * calculates the P3M @f$ \tilde{A}_2 @f$ term by
     * @f[ \tilde{A}_2 = \frac{4\pi}{V\tilde{S}^2} \sum_{\mathbf{l} \in G, l\neq 0}
             \hat{G}_g^{\dagger}(\mathbf{k}_{\mathbf{l}})
             \hat{R}_g(\mathbf{k}_{\mathbf{l}}) @f]
     *
     * it also calculates the derivative which is used to calculate the virial
     * @f[ \underline{\frac{d\tilde{A}_2}{dL}} = \frac{4\pi}{V\tilde{S}^2}
             \sum_{\mathbf{l} \in G, l\neq 0}
             \left[\hat{\underline{\Gamma}}_g^0(\mathbf{k}_{\mathbf{l}}) + 
             \mathrm{diag}(\hat{G}_g^{\dagger}(\mathbf{k}_{\mathbf{l}}))\right]
             \hat{R}_g(\mathbf{k}_{\mathbf{l}}) @f]
     *
     * @param[in] storage the storage is used to store the energy
     * @sa interaction::Nonbonded_Outerloop::ls_self_outerloop
     */
    template<class MeshType>
    static void calculate_p3m_selfterm(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim);
    
    /**
     * calculate the gridded electrostatic field
     *
     * From the fourier transform of the electrostatic potential the electrostatic
     * field is calculated by ik or finite differentiation. The electrostatic field is
     * stored in the configuration.
     *
     * The gridded electrostatic field is calculated either from the potential 
     * (which is itself calculated by a backward FFT) by finite differentiation
     * @f[ \mathbf{E}_g(\mathbf{r}_{\mathbf{n}}) =
               - \sum_{\mathbf{n}'\in G}\,i\,V_G\,
               \mathbf{D}_g(\mathbf{r}_{\mathbf{n}}-\mathbf{r}_{\mathbf{n}'})
               \Phi_{\gamma,g}\left(\mathbf{r}_{\mathbf{n}'}\right) @f]
     * or by multiplication by @f$ -i\mathbf{k}_{\mathbf{l}} @f$ by application of
     * @f[ \mathbf{\hat{E}}_g(\mathbf{r}_{\mathbf{n}}) = -i\mathbf{k}_{\mathbf{l}}
             \hat{\Phi}_{\gamma,g}\left(\mathbf{k}_{\mathbf{l}}\right) @f]
     * and application of three backards FFT to the components of the electrostatic
     * field.
     *
     * In the finite differentiation scheme the finite-difference operator
     * @f$ iV_G\mathbf{D}_g @f$ of order @f$ q @f$ estimates the gradient of
     * a given function at the @f$ 6q @f$ neighboring grid points along the box 
     * axes. 
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
     * The force is calculated by
     * @f[\mathbf{F}_{\gamma, i} = q_i E(\mathbf{r}_i) @f]
     * and
     * @f[\mathbf{E}(\mathbf{r}) = V_G \sum_{\mathbf{n} \in G} P(\mathbf{r} - 
              \mathbf{r}_\mathbf{n})\mathbf{E}_g(\mathbf{r}_{\mathbf{n}}) @f]
     *
     * In practice, summation is restricted over assignment-order neighboring
     * cells.
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
            std::vector<int> & domain,
            const math::VArray & r,
            const unsigned int size);
  };
 
  
}

#endif

