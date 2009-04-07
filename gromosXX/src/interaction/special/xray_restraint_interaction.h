/**
 * @file xray_restraint_interaction.h
 * xray restraining
 */

#ifndef INCLUDED_XRAY_RESTRAINT_INTERACTION_H
#define INCLUDED_XRAY_RESTRAINT_INTERACTION_H

// Additional Clipper Headers
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

namespace interaction {

  /**
   * @class xray_restraint_interaction
   * calculates the xray restraining interaction
   */
  class Xray_Restraint_Interaction : public Interaction {
  public:

    /**
     * Constructor.
     */
    Xray_Restraint_Interaction();
    /**
     * Destructor.
     */
    virtual ~Xray_Restraint_Interaction();

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

  protected:
    /**
     * pointer to the atoms
     */
    clipper::Atom_list atoms;
    clipper::HKL_info hkls;
    clipper::HKL_data<clipper::data32::F_phi> fphi;
    clipper::HKL_data<clipper::data32::F_phi> D_k;
    clipper::Xmap<clipper::ftype32> d_r;
    /**
     * decision-boolean for reseting averages
     */
    bool readavg;

    template<math::boundary_enum B, math::virial_enum V>
    void _calculate_xray_restraint_interactions
    (topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);



  };

} // interaction

#endif
