/**
 * @file molecular_virial_interaction.h
 * recover molecular virial from atomic virial
 */

#ifndef INCLUDED_MOLECULAR_VIRIAL_INTERACTION_H
#define INCLUDED_MOLECULAR_VIRIAL_INTERACTION_H

namespace interaction
{
  /**
   * @class molecular_virial_interaction
   * recovers a molecular virial from an atomic one
   */
  class Molecular_Virial_Interaction : public Interaction
  {
  public:
    /**
     * Constructor
     */
    Molecular_Virial_Interaction() : Interaction("MolecularVirial") {}
    
    /**
     * Destructor
     */
    virtual ~Molecular_Virial_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      if (!quiet)
	os << "Molecular virial (calculated from atomic virial)\n";
      return 0;
    };
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
  };
  
} // interaction

#endif
    
