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
      typedef MD<t_simulation, t_temperature, 
      t_pressure, t_distance_constraint,
      t_integration> parent_type;

      typedef interaction::Perturbed_Nonbonded_Inner_Loop<t_simulation,
      typename t_simulation::system_type> perturbed_innerloop_type;

      typedef interaction::Perturbed_Nonbonded_Inner_Loop_Virial<t_simulation,
      typename t_simulation::system_type> perturbed_innerloop_virial_type;
      
      typedef interaction::Basic_Pairlist<t_simulation,
      interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, 
      interaction::Twinrange_Chargegroup_Filter<t_simulation,
      interaction::Nonbonded_Base,
      interaction::Perturbed_Nonbonded_Inner_Loop<
      t_simulation, interaction::Storage>,
      interaction::Perturbation_Filter<t_simulation,
      interaction::Nonbonded_Base, true> > > >
      perturbed_pairlist_type;

      typedef interaction::Basic_Pairlist<t_simulation,
      interaction::Chargegroup_Range_Pairlist_Algorithm<t_simulation, 
      interaction::Twinrange_Chargegroup_Filter<t_simulation,
      interaction::Nonbonded_Base,
      interaction::Perturbed_Nonbonded_Inner_Loop_Virial<
      t_simulation, interaction::Storage>,
      interaction::Perturbation_Filter<t_simulation,
      interaction::Nonbonded_Base, true> > > > 
      perturbed_pairlist_virial_type;

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
				 io::InInput &input,
				 io::Argument &args);
    
      /**
       * initialize the perturbation parameters.
       */
      int init_perturbation(io::Argument &args, io::InInput &input);

      /**
       * print pairlists.
       */
      virtual void print_pairlists(); 
      
     /**
       * energy calculation.
       */
      virtual void do_energies();
    
    };
  
} // algorithm

#include "perturbation_md.tcc"

#endif
