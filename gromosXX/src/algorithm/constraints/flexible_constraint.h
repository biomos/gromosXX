/**
 * @file flexible_constraint.h
 * the flexible shake algorithm
 */

#ifndef INCLUDED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_FLEXIBLE_CONSTRAINT_H

namespace interaction
{
  class Nonbonded_Interaction;
}

namespace algorithm
{
  /**
   * @class Flexible_Constraint
   * calculates the flexible constraint distance
   */
  class Flexible_Constraint : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Flexible_Constraint(double const tolerance = 0.000001,
			int const max_iterations = 1000,
			interaction::Forcefield *ff = NULL);

    /**
     * Destructor.
     */
    virtual ~Flexible_Constraint();

    /**
     * initialization.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     bool quiet = false);

    /**
     * apply flexible shake.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    
    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * tolerance.
     */
    double const tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const max_iterations()const {return m_max_iterations;}

    /**
     * the const bond type parameter.
     */
    std::vector<interaction::bond_type_struct> const &parameter()const
    {
      return m_parameter;
    }
    /**
     * the bond type parameter.
     */
    std::vector<interaction::bond_type_struct> & parameter()
    {
      return m_parameter;
    }

    void calc_distance
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim
     );

    void solute
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     int & error
     );

    void _store_lengths
    (
     configuration::Configuration & conf
     );

  protected:
    /**
     * shake tolerance
     */
    double m_tolerance;
    /**
     * max iterations.
     */
    const int m_max_iterations;
    /**
     * bond parameter
     */
    std::vector<interaction::bond_type_struct> m_parameter;

    /**
     * flexible constraint lengths
     */
    std::vector<double> m_flex_len;
    
    /**
     * the nonbonded's for exact algorithm...
     * (only implemented for diatomic molecules...
     *  no bonded terms!)
     */
    interaction::Nonbonded_Interaction * m_nonbonded;

    /**
     * the (undetermined) constrained forces (without the lambda factor)
     */
    std::vector<std::vector<math::Vec> > m_force;

    template<math::boundary_enum B, math::virial_enum V>
    int _iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     bool & convergence,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     double dt,
     math::Periodicity<B> const & periodicity,
     bool do_constraint_force = false
     );

    template<math::boundary_enum B, math::virial_enum V>
    int _exact_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     bool & convergence,
     double dt,
     math::Periodicity<B> const & periodicity
     );

    template<math::boundary_enum B, math::virial_enum V>
    void _calc_distance
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim
     );

    template<math::boundary_enum B, math::virial_enum V>
    void _calc_undetermined_forces
    (
     topology::Topology &topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

    template<math::boundary_enum B, math::virial_enum V>
    void _solute
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     int & error
     );

    template<math::boundary_enum B, math::virial_enum V>
    int _approx_work
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     double dt
     );

  };
  
} // algorithm

#endif
