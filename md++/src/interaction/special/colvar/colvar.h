/**
 * @file colvar.h
 * @brief Base class for collective variables (colvars)
 */

#ifndef INCLUDED_COLVAR_H
#define INCLUDED_COLVAR_H

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
    Colvar(std::string name) : name(name), w0(1), w0A(1), w0B(1), m_timer(name) {};
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
    /**
     * store time used in algorithm.
     */
    util::Algorithm_Timer m_timer;
    
  };  
  
} // interaction

#endif
