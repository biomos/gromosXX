/**
 * @file settle.h
 * the settle algorithm.
 */

#ifndef INCLUDED_SETTLE_H
#define INCLUDED_SETTLE_H

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm {

  /**
   * @class Settle
   * implements the settle algorithm for 3 site water models.
   *
   * Coordinates and velocities are reset according to
   * S. Miyamoto and P. A. Kollman, SETTLE: An Analytical Version of the SHAKE
   * and RATTLE Algorithm for Rigid Water Models, J. Comput. Chem 13, Issue 8, 
   * 1992, pp. 952-962
   *
   * The source code is annotaed according to the variables names used in the 
   * paper. 
   *
   * The constraint force and virial is not calculated as discribed in the paper
   * but from the displacement due to the constraint.
   */
  class Settle : public Algorithm {
  public:

    /**
     * Constructor.
     */
    Settle() : Algorithm("Settle") {
    }

    /**
     * Destructor.
     */
    virtual ~Settle() {
    }

    /**
     * apply the SETTLE algorithm
     */
    virtual int apply(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);

    /**
     * initialize startup positions and velocities
     * if required.
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            std::ostream & os = std::cout,
            bool quiet = false);

    /**
     * the const bond type parameter.
     */
    std::vector<interaction::bond_type_struct> const &parameter()const {
      return m_parameter;
    }

    /**
     * the bond type parameter.
     */
    std::vector<interaction::bond_type_struct> & parameter() {
      return m_parameter;
    }

    /**
     * accessor to the constrained atoms
     */
    std::set<unsigned int> & constrained_atoms() {
      return m_constrained_atoms;
    }

    /**
     * accessor to the constrained atoms
     */
    const std::set<unsigned int> & constrained_atoms() const {
      return m_constrained_atoms;
    }

  protected:
    /**
     * print an error message to std::cout if SETTLE fails
     */
    void printError(topology::Topology const & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            unsigned int atom, std::string message);

    /** 
     * apply it for the solvent with correct periodic boundary conditions
     */
    void solvent(topology::Topology const & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            int & error);

    /**
     * bond parameter
     */
    std::vector<interaction::bond_type_struct> m_parameter;
    /**
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;
    /** 
     * rank and size for parallelization
     */
    int m_rank, m_size;
  };

} //algorithm

#endif
