/**
 * @file contactnum.h
 * contact number restraining
 */

#ifndef INCLUDED_CONTACTNUM_RESTRAINT_INTERACTION_H
#define INCLUDED_CONTACTNUM_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class contactnum_colvar
   * calculates contact number and derivatives with respect to the position
   * which can then be used for energy and force calculation in the colvar
   * restraint interaction
   */
  class Contactnum_Colvar : public Colvar
  {
  public:
    /**
     * Constructor.
     */
    Contactnum_Colvar() : Colvar("Contacnum") {}
    
    /**
     * Destructor.
     */
    virtual ~Contactnum_Colvar() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    /**
     * calculate contact number and derivatives
     */
    virtual int calculate(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
				       
    topology::contactnum_restraint_struct *params;
    
    private:
      int mm, nn;
      double rcut;    
      double switchingfunction(double rdist,double&dfunc,int nn,int mm);
      double fastpow(double base, int exp);
    
  };
  
} // interaction

#endif
