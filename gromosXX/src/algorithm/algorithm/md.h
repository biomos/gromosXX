/**
 * @file md.h
 * the md algorithm.
 */

#ifndef INCLUDED_MD_H
#define INCLUDED_MD_H

namespace algorithm
{
  /**
   * @class MD
   * the MD algorithm
   */
  template<typename t_md_spec=MD_spec<>,
	   typename t_interaction_spec=Interaction_spec<
    typename t_md_spec::simulation_type>
  >
  class MD : public MD_Base<t_md_spec, t_interaction_spec>
  {
  public:
    /**
     * Constructor.
     */
    MD();

    /**
     * Destructor.
     */
    virtual ~MD();
    
    /**
     * perform an md simulation.
     * calls run.
     */
    virtual int do_md(io::Argument &args, io::InInput &input);

    /**
     * run the system.
     * @param time the time to run the system.
     */
    virtual void run(double time = -1);
    
    /**
     * create a Gromos96 like forcefield.
     */
    virtual void G96Forcefield(io::InTopology &topo,
			       io::InInput &input,
			       io::Argument &args);
      
  protected:
    
    /**
     * open the input files.
     */
    virtual void open_files(io::Argument &args, io::InTopology &topo,
		    io::InTrajectory &sys, io::InInput &input);
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
    /**
     * initialize the output.
     */
    virtual void init_output(io::Argument &args, io::InInput &input);

    /**
     * initialization.
     */
    virtual int initialize(io::Argument &args, io::InInput &input);
    
    /**
     * pre md
     */
    virtual int pre_md(io::InInput &input);

    /**
     * before we do a step in run md
     */
    virtual void pre_step();
    
    /**
     * do the step in run md.
     */
    virtual void do_step();
    
    /**
     * after the step in run md.
     */
    virtual void post_step();

    /**
     * after the md run.
     */
    virtual void post_md();

    /**
     * print pairlists.
     */
    virtual void print_pairlists();
    /**
     * calculate and print the energies.
     */
    virtual void do_energies();
    
  };
  
} // algorithm

#include "md.tcc"

#endif
