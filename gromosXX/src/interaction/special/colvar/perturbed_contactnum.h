/**
 * @file perturbed_contactnum.h
 * perturbed contact number restraining
 */

#ifndef INCLUDED_PERTCONTACTNUM_RESTRAINT_INTERACTION_H
#define INCLUDED_PERTCONTACTNUM_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Perturbed_Contactnum_Colvar
   * calculates contact number and contact number derivatives with respect to the position
   * which can then be used for energy and force calculation in the colvar
   * restraint interaction
   */
  class Perturbed_Contactnum_Colvar : public Colvar
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Contactnum_Colvar() : Colvar("PerturbedContactnum") {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Contactnum_Colvar() {}

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
				       
    topology::perturbed_contactnum_restraint_struct *params;
    
    private:
      int mm, nn;
      double rcut;    
      double switchingfunction(double rdist,double&dfunc,int nn,int mm);
      double fastpow(double base, int exp);
    
  };
  
} // interaction

#endif
