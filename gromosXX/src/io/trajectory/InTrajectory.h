/**
 * @file InTrajectory.h
 * read in a G96 trajectory file.
 */

#ifndef INCLUDED_INTRAJECTORY_H
#define INCLUDED_INTRAJECTORY_H

namespace io {

  /**
   * @class InTrajectory
   * reads in a trajectory file and parses
   * it into simulation::system.
   */
  class InTrajectory : public GInStream {

  public:
    /**
     * Constructor.
     */
    InTrajectory(std::istream& is);
    /**
     * Read in a G96 trajectory into the system.
     */
    InTrajectory & operator>>(simulation::system &sys);
    /**
     * switch: read in position (positionred) or not.
     */
    bool read_position;
    /**
     * switch: read in velocity (velocityred) or not.
     */
    bool read_velocity;
    /**
     * switch: read in box or not.
     */
    bool read_box;

  private:
    /**
     * read POSITION block.
     */
    bool _read_position(math::VArray &pos, std::vector<std::string> &buffer);
    /**
     * read POSITIONRED block.
     */
    bool _read_positionred(math::VArray &pos, std::vector<std::string> &buffer);
    /**
     * read VELOCITY block.
     */
    bool _read_velocity(math::VArray &vel, std::vector<std::string> &buffer);
    /**
     * read POSITIONRED block.
     */
    bool _read_velocityred(math::VArray &vel, std::vector<std::string> &buffer);
    /**
     * read BOX block.
     */
    bool _read_box(simulation::system &sys, std::vector<std::string> &buffer);
    
  };
  

} // io

// template and inline methods
#include "InTrajectory.tcc"

#endif
