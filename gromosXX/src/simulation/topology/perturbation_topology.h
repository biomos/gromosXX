/**
 * @file perturbation_topology.h
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
     * perturbed solute accessor as const.
     */
    Perturbed_Solute const & perturbed_solute()const;
    
    /**
     * lambda accessor.
     */
    double lambda();
    /**
     * lambda accessor const.
     */
    double const lambda()const;
    
    /**
     * lambda accessor.
     */
    void lambda(double const l);

    /**
     * alpha for lennard jones
     */
    double alpha_lj();
    
    /**
     * alpha for lennard jones
     */
    void alpha_lj(double const a);
    
    /**
     * alpha for coulomb
     */
    double alpha_crf();
    /**
     * alpha for coulomb
     */
    void alpha_crf(double const a);
    
    /**
     * lambda dependence for nonbonded interactions
     */
    int nlam();
    /**
     * lambda dependence for nonbonded interactions
     */
    void nlam(int const n);
    
    /**
     * const perturbed atom accessor. (bool)
     */
    std::vector<bool> const & perturbed_atom()const;

    /**
     * perturbed atom accessor. (bool)
     */
    std::vector<bool> & perturbed_atom();

    /**
     * recalculate lambda dependent properties.
     */
    void update_for_lambda();

    /**
     * calculate constraint degrees of freedom
     * also adding the perturbed constraints!
     */
    void calculate_constraint_dof(simulation::Multibath &multibath)const;
    
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
    /**
     * alpha lennard jones
     */
    double m_alpha_lj;
    /**
     * alpha coulomb
     */
    double m_alpha_crf;
    /**
     * lambda dependence
     */
    int m_nlam;
    
    
  };

} // simulation

// inline methods6

#include "perturbation_topology.tcc"

#endif
