/**
 * @file perturbation_md.h
 * md algorithm including perturbation.
 * should maybe be changed to a mix in class
 * using multiple inheritance???
 */

#ifndef INCLUDED_PERTURBATION_MD_H
#define INCLUDED_PERTURBATION_MD_H

namespace algorithm
{
  /**
   * @class Perturbation_MD
   * MD algorithm with perturbation.
   */
  template<typename t_simulation,
	   typename t_temperature = algorithm::Berendsen_Thermostat,
	   typename t_pressure = algorithm::Berendsen_Barostat,
	   typename t_distance_constraint = algorithm::Shake<t_simulation>,
	   typename t_integration = algorithm::Leap_Frog<t_simulation> >
  class Perturbation_MD : public MD<t_simulation, t_temperature,
				    t_pressure, t_distance_constraint,
				    t_integration>
  {
  public:
    /**
     * Constructor.
     */
    Perturbation_MD(t_simulation &sim);

  protected:
    /**
     * initialize the input.
     */
    virtual void init_input(io::Argument &args, io::InTopology &topo,
			    io::InTrajectory &sys, io::InInput &input);

    /**
     * read the input and setup a standard simulation.
     */
    virtual void read_input(io::Argument &args, io::InTopology &topo,
			    io::InTrajectory &sys, io::InInput &input);
    
    virtual void G96Forcefield(io::InTopology &topo,
			       io::InInput &input);
    
    /**
     * initialize the perturbation parameters.
     */
    int init_perturbation(io::Argument &args, io::InInput &input);
    
  };
  
} // algorithm

#include "perturbation_md.tcc"

#endif
