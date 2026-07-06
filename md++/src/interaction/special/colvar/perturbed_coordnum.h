/**
 * @file perturbed_coordnum.h
 * @brief Perturbed coordination number colvar
 */

#ifndef INCLUDED_PERTURBED_COORDNUM_H
#define INCLUDED_PERTURBED_COORDNUM_H

namespace interaction {

  class Perturbed_Coordnum_Colvar : public Colvar {
  public:
    Perturbed_Coordnum_Colvar() : Colvar("PerturbedCoordnum") {}
    virtual ~Perturbed_Coordnum_Colvar() {}

    virtual int init(topology::Topology &topo,
                     configuration::Configuration &conf,
                     simulation::Simulation &sim,
                     std::ostream &os = std::cout,
                     bool quiet = false);

    virtual int calculate(topology::Topology &topo,
                          configuration::Configuration &conf,
                          simulation::Simulation &sim);

    topology::perturbed_coordnum_restraint_struct *params;

  private:
    int mm, nn;
    double rcut;
  };

} // namespace interaction

#endif
