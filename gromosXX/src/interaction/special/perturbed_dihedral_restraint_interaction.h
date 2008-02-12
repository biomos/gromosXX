/**
 * @file perturbed_dihedral_restraint_interaction.h
 * perturbed dihedral restraining
 */

#ifndef INCLUDED_PERTURBED_DIHEDRAL_RESTRAINT_INTERACTION_H
#define INCLUDED_PERTURBED_DIHEDRAL_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Perturbed_Dihedral_Restraint_Interaction
   * calculates the perturbed dihedral restraining interaction
   */
  class Perturbed_Dihedral_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Dihedral_Restraint_Interaction() : Interaction("PerturbedDihedralRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Dihedral_Restraint_Interaction() {}

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
	os << "Perturbed dihedral restraint interaction\n";
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
