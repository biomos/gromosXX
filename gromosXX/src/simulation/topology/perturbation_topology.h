/**
 * @file perturbation_topology_h
 * a perturbation topology
 * (also contains the unperturbed topology)
 */

#ifndef INCLUDED_PERTURBATION_TOPOLOGY_H
#define INCLUDED_PERTURBATION_TOPOLOGY_H

namespace simulation
{

  /**
   * @class Perturbation_Topology
   * holds the topological information
   * of the end states of a perturbed 
   * simulation. (plus the normal
   * topology)
   */
  class Perturbation_Topology : public Topology
  {
  public:
    /**
     * Constructor.
     */
    explicit Perturbation_Topology();

    /**
     * set the capacity of the arrays.
     * @override
     */
    void resize(size_t const atoms);

    /**
     * perturbed solute accessor.
     */
    Perturbed_Solute &perturbed_solute();

    /**
     * lambda accessor.
     */
    double lambda();

    /**
     * lambda accessor.
     */
    void lambda(double const l);

    /**
     * const perturbed atom accessor. (bool)
     */
    std::vector<bool> const & perturbed_atom()const;

    /**
     * perturbed atom accessor. (bool)
     */
    std::vector<bool> & perturbed_atom();
    
  private:
    /**
     * is the atom perturbed?
     */
    std::vector<bool> m_perturbed_atom;
    
    /**
     * the perturbed solute information.
     */
    Perturbed_Solute m_perturbed_solute;

    /**
     * lambda.
     */
    double m_lambda;
    
  };

} // simulation

// inline methods6

#include "perturbation_topology.tcc"

#endif
