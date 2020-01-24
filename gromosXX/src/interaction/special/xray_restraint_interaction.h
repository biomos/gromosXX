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
            topology::xray_umbrella_weight_struct & param,
            configuration::Configuration & conf,
            clipper::Atom_list & atoms,
            clipper::Xmap<clipper::ftype32> & rho_calc,
            clipper::Xmap<clipper::ftype32> & rho_obs,
            double to_ang) :
            weight(0.0), param(param),
                    conf(conf), atoms(atoms), rho_calc(rho_calc), rho_obs(rho_obs), to_ang(to_ang) {
    }
#endif

    virtual double get_weight() const { return weight; }
    virtual void increment_weight();
    virtual void write(std::ostream & os) const;
    virtual void read(std::istream & is) { is >> weight; }
  protected:
    double weight;
#ifdef HAVE_CLIPPER
    topology::xray_umbrella_weight_struct & param;
    configuration::Configuration & conf;
    clipper::Atom_list & atoms;
    clipper::Xmap<clipper::ftype32> & rho_calc;
    clipper::Xmap<clipper::ftype32> & rho_obs;
    double to_ang;
#endif
  };

/**
   * electron density umbrella weight factory
   */
  class Electron_Density_Umbrella_Weight_Factory : public util::Umbrella_Weight_Factory {
  public:
#ifdef HAVE_CLIPPER
    Electron_Density_Umbrella_Weight_Factory(
            topology::xray_umbrella_weight_struct & param,
            configuration::Configuration & conf,
            clipper::Atom_list & atoms,
            clipper::Xmap<clipper::ftype32> & rho_calc,
            clipper::Xmap<clipper::ftype32> & rho_obs,
            double to_ang) :
    param(param), conf(conf),
    atoms(atoms), rho_calc(rho_calc), rho_obs(rho_obs), to_ang(to_ang) {}
#endif
    virtual util::Umbrella_Weight * get_instance();
#ifdef HAVE_CLIPPER
  private:
    topology::xray_umbrella_weight_struct & param;
    configuration::Configuration & conf;
    clipper::Atom_list & atoms;
    clipper::Xmap<clipper::ftype32> & rho_calc;
    clipper::Xmap<clipper::ftype32> & rho_obs;
    double to_ang;
 #endif
  };

} // interaction
#endif

