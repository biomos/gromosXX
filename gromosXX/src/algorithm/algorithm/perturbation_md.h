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
  template<typename t_spec = perturbed_MD_spec>
    class Perturbation_MD : public MD<t_spec>
    {
      public:
      /**
       * Constructor.
       */
      Perturbation_MD();

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
