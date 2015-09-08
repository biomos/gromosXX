/**
 * @file rdc_restraint_interaction.h
 * RDC restraining
 */

#ifndef __RDC_RESTRAINT_INTERACTION_H__
#define __RDC_RESTRAINT_INTERACTION_H__

#include <stdheader.h>
#include <topology/topology.h>
#include <interaction/interaction.h>

//  //math::eps0_i is defined as N_A/\epsilon_0 = 1.7459 *10^3 (kJ nm)/(e^2 mol)
//  //math::spd_l is defined as c = 2.9979 *10^5 nm/ps
//  const double m0 = math::eps0_i / (math::spd_l * math::spd_l); //1.9426 *10^-8 (kJ ps^2)/(e^2 mol nm)
//  //math::h_bar is defined as  hbar*N_A = 6.351 *10^-2 (kJ ps)/mol
//  const double h_planck = math::h_bar * 2.0 * math::Pi; // 0.399 (kJ ps)/mol


namespace interaction {

  // returns rdc_max * r^3 (because 1/r^3 needs to be averaged) in units of nm^3/ps
  inline double RDC_max(std::vector<topology::rdc_restraint_struct>::iterator it){
    return - (math::eps0_i * math::h_bar * it->gyri * it->gyrj) / pow(2.0 * math::Pi * math::spd_l, 2);
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
     * Default Constructor
     */
    RDC_Restraint_Interaction(): Interaction("RDCRestraint") {}

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

  };
} // namespace interaction

#endif	// __RDC_RESTRAINT_INTERACTION_H__

