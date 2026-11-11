/**
 * @file colvar.h
 * @brief Base class for collective variables (colvars)
 */

#ifndef INCLUDED_COLVAR_H
#define INCLUDED_COLVAR_H

#include "../../../stdheader.h"

namespace configuration{
  class Configuration;
}
namespace topology{
  class Topology;
}
namespace simulation{
  class Simulation;
}
namespace util {
  class Algorithm_Timer;
  class Virtual_Atom;
}

namespace interaction {
  /**
   * @class Colvar
   * @interface Colvar
   * declares the colvar interface.
   */
  class Colvar
  {
  public:
    /**
     * Constructor.
     */
    Colvar(std::string name) : name(name), w0(1), w0A(1), w0B(1), strain_derivative(0.0), m_timer(name) {};
    /**
     * Destructor.
     */
    virtual ~Colvar(){};
    /**
     * the name of the colvar.
     * can be used to identify a special class.
     */
    std::string name;
    
    /**
     * current value
     */
    double ct;
    /**
     * target value and force constant weight
     */
    double targetvalue, w0;
    /**
     * A and B state target value and force constant weight for the perturbed colvars
     */
    double targetvalueA, targetvalueB, w0A, w0B;
    /**
     * pointers to all atoms for which we need to get forces
     */
    std::vector<util::Virtual_Atom* > atoms;
    /**
     * derivatives with respect to the atom positions for the above atoms
     */
    math::VArray derivatives;

    /**
     * Derivative of the CV under affine deformation of its endpoint positions
     * and periodic box: component (a,b) corresponds to x_b += epsilon * x_a.
     * Geometry must use the same periodic images as the value calculation.
     */
    math::Matrix strain_derivative;

    /**
     * Apply the bias forces through the endpoint coordinate derivatives.
     * Optionally add the geometry virial and virtual-site chain-rule correction.
     * molecular_virial converts this contribution to the native pressure-group
     * convention without including pre-existing forces.
     */
    void apply_forces(topology::Topology &topo, configuration::Configuration &conf,
                      double dE_dinstant, bool include_virial,
                      bool molecular_virial=false) const;
    
    /** Apply into state-specific buffers, e.g. before EDS mixing. */
    void apply_forces(topology::Topology &topo, configuration::Configuration &conf,
                      double dE_dinstant, bool include_virial,
                      math::VArray &forces, math::Matrix &virial,
                      bool molecular_virial=false) const;

    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) = 0;
    
    /**
     * calculate the collective variable and derivatives.
     */
    virtual int calculate(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim) = 0;

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      m_timer.print(os);
    }
    
  protected:
    /** Add the affine-endpoint part of the internal virial. */
    void add_virial(math::Matrix &virial, double dE_dinstant) const;

    /**
     * store time used in algorithm.
     */
    util::Algorithm_Timer m_timer;
    
  };  
  
} // interaction

#endif
