/**
 * @file dihedral_restraint_interaction.h
 * dihedral restraining
 */

#ifndef INCLUDED_DIHEDRAL_RESTRAINT_INTERACTION_H
#define INCLUDED_DIHEDRAL_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class dihedral_restraint_interaction
   * calculates the dihedral restraining interaction
   */
  class Dihedral_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Dihedral_Restraint_Interaction() : Interaction("DihedralRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Dihedral_Restraint_Interaction() {}

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
	os << "Dihedral restraint interaction\n";
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
