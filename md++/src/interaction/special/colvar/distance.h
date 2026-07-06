/**
 * @file distance.h
 * @brief Distance collective variable
 */

#ifndef INCLUDED_DISTANCE_COLVAR_H
#define INCLUDED_DISTANCE_COLVAR_H

namespace interaction {

  /**
   * @class Distance_Colvar
   * @brief Calculates a distance collective variable and derivatives.
   *
   * RAH dimensionality is handled here because it changes the definition of
   * the distance itself.  RAH half-harmonic behaviour is exposed via
   * rah_mode_value() and handled by the generic Colvar_Bias layer.
   */
  class Distance_Colvar : public Colvar {
  public:
    Distance_Colvar() : Colvar("Distance"), rah_mode(0) {}
    virtual ~Distance_Colvar() {}

    virtual int init(topology::Topology &topo,
                     configuration::Configuration &conf,
                     simulation::Simulation &sim,
                     std::ostream &os = std::cout,
                     bool quiet = false);

    virtual int calculate(topology::Topology &topo,
                          configuration::Configuration &conf,
                          simulation::Simulation &sim);

    int rah_mode_value() const { return rah_mode; }

    topology::distance_restraint_struct *params;

  private:
    math::Vec dim_mask;
    int rah_mode;
  };

} // namespace interaction

#endif