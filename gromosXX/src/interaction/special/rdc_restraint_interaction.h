//GP AT WORK 23
/**
 * @file rdc_restraint_interaction.h
 * RDC restraining
 */

#ifndef __RDC_RESTRAINT_INTERACTION_H__
#define __RDC_RESTRAINT_INTERACTION_H__

/* Define the way the atoms should be represented.
  - Magnetic field vectors
  - Alignment tensor with 5 components
  - Spherical harmonics
*/

#include "../../math/random.h"

//  //math::eps0_i is defined as N_A/\epsilon_0 = 1.7459 *10^3 (kJ nm)/(e^2 mol)
//  //math::spd_l is defined as c = 2.9979 *10^5 nm/ps
//  const double m0 = math::eps0_i / (math::spd_l * math::spd_l); //1.9426 *10^-8 (kJ ps^2)/(e^2 mol nm)
//  //math::h_bar is defined as  hbar*N_A = 6.351 *10^-2 (kJ ps)/mol
//  const double h_planck = math::h_bar * 2.0 * math::Pi; // 0.399 (kJ ps)/mol


namespace interaction {

  // returns rdc_max * r^3 (because 1/r^3 needs to be averaged) in units of nm^3/ps
  inline double RDC_max(std::vector<topology::rdc_restraint_struct>::iterator it){
    //return - (math::eps0_i * math::h_bar * it->gyri * it->gyrj) / (pow(math::spd_l,2) * 4.0 * pow(math::Pi,2));
    return - (math::eps0_i * math::h_bar * it->gyri * it->gyrj) / pow(2.0 * math::Pi * math::spd_l,2);
    // keeping in mind: ( kJ^2 * ps^3 )/( u^2 * mol^2 * nm ) * 1/N_A^2 = 1 nm^3/ps
    // i.e. we implicitly divide by N_A^2 but the factor is absorbed in units ...
  }

  /**
   * @class RDC_Restraint_Interaction
   * calculates the RDC restraining interaction
   */
  class RDC_Restraint_Interaction: public Interaction{
  private:

  public:
    /**
     * Constructor.
     */
    RDC_Restraint_Interaction(): Interaction("RDCRestraint") {}

    /**
     * Destructor.
     */
    virtual ~RDC_Restraint_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);


    /*
    template<math::boundary_enum B, math::virial_enum V> int
    _calculate_interactions_mfield(topology::Topology & topo,
                                   configuration::Configuration & conf,
                                   simulation::Simulation & sim);

    template<math::boundary_enum B, math::virial_enum V> int
    _calculate_interactions_tensor(topology::Topology & topo,
                                   configuration::Configuration & conf,
                                   simulation::Simulation & sim);

    template<math::boundary_enum B, math::virial_enum V> int
    _calculate_interactions_spherical(topology::Topology & topo,
                                   configuration::Configuration & conf,
                                   simulation::Simulation & sim);


    virtual double _calculate_derivative
    (topology::Topology & topo,
                                         configuration::Configuration &conf,
                                         simulation::Parameter const & param,
                                         std::vector<topology::rdc_restraint_struct>::const_iterator it,
                                         unsigned int eps_i,
                                         double RDCcurr, double RDCav,
                                         double theta);
    */


    /*
    struct sd_struct {
        math::SArray gammas;
        math::SArray c1, c2, c3, c4, c5, c6, c7, c8, c9;
        **
         * resize the variables
         *
        void resize(int size) {
            gammas.resize(size);
            c1.resize(size);
            c2.resize(size);
            c3.resize(size);
            c4.resize(size);
            c5.resize(size);
            c6.resize(size);
            c7.resize(size);
            c8.resize(size);
            c9.resize(size);
        }

        **
         * clear the variables
         *
        void clear() {
            gammas.clear();
            c1.clear();
            c2.clear();
            c3.clear();
            c4.clear();
            c5.clear();
            c6.clear();
            c7.clear();
            c8.clear();
            c9.clear();
        }
        **
         * constructor
         *
        sd_struct() :
            gammas(0), c1(0), c2(0), c3(0), c4(0), c5(0), c6(0), c7(0), c8(0), c9(0) {}
    };*/
  };
} // namespace interaction

#endif	// __RDC_RESTRAINT_INTERACTION_H__

