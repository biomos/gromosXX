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
     * read step block.
     */
    void read_STEP(int &num_steps, double &t0, double &dt);

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
