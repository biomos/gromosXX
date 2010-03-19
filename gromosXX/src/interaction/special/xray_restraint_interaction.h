/**
 * @file xray_restraint_interaction.h
 * xray restraining
 */

#ifndef INCLUDED_XRAY_RESTRAINT_INTERACTION_H
#define INCLUDED_XRAY_RESTRAINT_INTERACTION_H

// Additional Clipper Headers
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#endif

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
#ifdef HAVE_CLIPPER
    /**
     * the atoms
     */
    clipper::Atom_list atoms;
    /**
     * the calculated electron density
     */
    clipper::Xmap<clipper::ftype32> rho_calc;
    /**
     * the observed electron density
     */
    clipper::Xmap<clipper::ftype32> rho_obs;
    /**
     * the HKLs (reflections)
     */
    clipper::HKL_info hkls;
    /**
     * the structure factors
     */
    clipper::HKL_data<clipper::data32::F_phi> fphi;
    /**
     * structure factos built from the observed amplitudes and the
     * calculated phases
     */
    clipper::HKL_data<clipper::data32::F_phi> fphi_obs;
    /**
     * the gradient map
     */
    clipper::FFTmap_p1 D_k;
    /**
     * the map for the gradient convolution
     */
    clipper::Xmap<clipper::ftype32> d_r;
    /**
     * spacegroup for NCS restraints
     */
    clipper::Spacegroup ncs_spacegroup;
#endif

    template<math::boundary_enum B, math::virial_enum V>
    void _calculate_xray_restraint_interactions
    (topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            int & error);
  };

  /**
   * electron density umbrella weight
   */
  class Electron_Density_Umbrella_Weight : public util::Umbrella_Weight {
  public:
#ifdef HAVE_CLIPPER
    Electron_Density_Umbrella_Weight(
            std::vector<unsigned int> & variable_atoms,
            double threshold, double cutoff,
            configuration::Configuration & conf,
            clipper::Atom_list & atoms,
            clipper::Xmap<clipper::ftype32> & rho_calc,
            clipper::Xmap<clipper::ftype32> & rho_obs,
            double to_ang) :
            weight(0.0), variable_atoms(variable_atoms), threshold(threshold), cutoff(cutoff),
                    conf(conf), atoms(atoms), rho_calc(rho_calc), rho_obs(rho_obs), to_ang(to_ang) {
    }
#endif
    virtual double get_weight() const { return weight; }
    virtual void increment_weight();
    virtual void write(std::ostream & os) const;
    virtual void read(std::istream & is) { is >> weight; }
  protected:
    double weight;
    std::vector<unsigned int> & variable_atoms;
    double threshold;
    double cutoff;
    configuration::Configuration & conf;

#ifdef HAVE_CLIPPER
    clipper::Atom_list & atoms;
    clipper::Xmap<clipper::ftype32> & rho_calc;
    clipper::Xmap<clipper::ftype32> & rho_obs;
#endif
    double to_ang;
  };

/**
   * electron density umbrella weight factory
   */
  class Electron_Density_Umbrella_Weight_Factory : public util::Umbrella_Weight_Factory {
  public:
#ifdef HAVE_CLIPPER
    Electron_Density_Umbrella_Weight_Factory(std::vector<unsigned int> variable_atoms,
            double threshold, double cutoff,
            configuration::Configuration & conf,
            clipper::Atom_list & atoms,
            clipper::Xmap<clipper::ftype32> & rho_calc,
            clipper::Xmap<clipper::ftype32> & rho_obs,
            double to_ang) :
    variable_atoms(variable_atoms), threshold(threshold), cutoff(cutoff), conf(conf),
    atoms(atoms), rho_calc(rho_calc), rho_obs(rho_obs), to_ang(to_ang) {
    }
#endif
    virtual util::Umbrella_Weight * get_instance();
  protected:
    std::vector<unsigned int> variable_atoms;
    double threshold;
    double cutoff;
    configuration::Configuration & conf;
#ifdef HAVE_CLIPPER
    clipper::Atom_list & atoms;
    clipper::Xmap<clipper::ftype32> & rho_calc;
    clipper::Xmap<clipper::ftype32> & rho_obs;
#endif
    double to_ang;
  };

} // interaction
#endif

