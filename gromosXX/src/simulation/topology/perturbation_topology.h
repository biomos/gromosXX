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

  private:
    /**
     * is the atom perturbed?
     */
    std::vector<bool> m_perturbed_atom;
  };

} // simulation

// inline methods6

#include "perturbation_topology.tcc"

#endif
