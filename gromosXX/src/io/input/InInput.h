/**
 * @file InInput.h
 * read in a G96 input file.
 */

#ifndef INCLUDED_ININPUT_H
#define INCLUDED_ININPUT_H

namespace io {

  /**
   * @class InInput
   * reads in an input file and parses it.
   */
  class InInput : public GInStream {

  public:
    /**
     * Constructor.
     * read in the complete file at once.
     */
    InInput(std::istream& is) : GInStream(is) { read_stream(); };

    /**
     * Store standard parameters in the simulation.
     */
    template<typename t_topology, typename t_system>
    InInput & operator>>(simulation::Simulation<t_topology, t_system> &sim);
    
    /**
     * read step block.
     */
    void read_STEP(int &num_steps, double &t0, double &dt);

    /**
     * read the SYSTEM block.
     */
    void read_SYSTEM(int &nsm);
    
    /**
     * read the SHAKE block.
     */
    void read_SHAKE(int &ntc, double &tolerance);
    /**
     * read FORCE block.
     */
    void read_FORCE(int &do_bond, int &do_angle, int &do_improper,
		    int &do_dihedral, int &do_nonbonded);
    /**
     * read PRINT and WRITE block.
     */
    void read_PRINT(int &print_trajectory, int &print_velocity,
		    int &print_energy);

  private:
    /**
     * read the entire stream and store the blocks in the map.
     */
    void read_stream();

    std::map<std::string, std::vector<std::string> > m_block;
    
  };
  
} // io

// template and inline methods
#include "InInput.tcc"

#endif
