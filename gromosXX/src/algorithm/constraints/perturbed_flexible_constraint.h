/**
 * @file perturbed_flexible_constraint.h
 * the perturbed flexible constraint algorithm.
 */

#ifndef INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H

namespace algorithm
{
  /**
   * @class Perturbed_Flexible_Constraint
   * implements the flexible constraint algorithm 
   * for perturbed distance constraints.
   */
  class Perturbed_Flexible_Constraint : public Flexible_Constraint
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Flexible_Constraint
    (
     double const tolerance = 0.000001,
     int const max_iterations = 1000,
     interaction::Forcefield *ff = NULL
     );
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Flexible_Constraint();
    
    /**
     * apply perturbed flexible constraints.
     * also applies unperturbed constraints
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    
    /**
     * initialize startup positions and velocities
     * if required.
     */
    int init(topology::Topology & topo,
	     configuration::Configuration & conf,
	     simulation::Simulation & sim,
	     bool quiet = false);

  protected:
  private:
    template<math::boundary_enum B, math::virial_enum V>
    int iteration
    (
     topology::Topology &topo,
     configuration::Configuration & conf,
     bool & convergence,
     std::vector<double> & flex_len,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     double dt,
     math::Periodicity<B> const & periodicity,
     bool do_constraint_force = false
     );

    template<math::boundary_enum B, math::virial_enum V>
    void calc_distance
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim,
     std::vector<double> & flex_len
     );

    template<math::boundary_enum B, math::virial_enum V>
    void solute
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     int & error
     );

  };
  
} //algorithm

#endif
