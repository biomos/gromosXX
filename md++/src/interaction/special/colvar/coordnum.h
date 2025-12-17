/**
 * @file coordnum.h
 * @brief Coordination number collective variable
 */

#ifndef INCLUDED_COORDNUM_RESTRAINT_INTERACTION_H
#define INCLUDED_COORDNUM_RESTRAINT_INTERACTION_H

namespace interaction {

  /**
   * @class Coordnum_Colvar
   * @brief Calculates coordination number and derivatives with respect to position.
   */
  class Coordnum_Colvar : public Colvar {
  public:
    Coordnum_Colvar() : Colvar("Coordnum") {}
    virtual ~Coordnum_Colvar() {}

    virtual int init(topology::Topology &topo,
                     configuration::Configuration &conf,
                     simulation::Simulation &sim,
                     std::ostream &os = std::cout,
                     bool quiet = false);

    virtual int calculate(topology::Topology &topo,
                          configuration::Configuration &conf,
                          simulation::Simulation &sim);

    topology::coordnum_restraint_struct *params;

  private:
    int mm, nn;
    double rcut;
  };

} // namespace interaction

#endif
