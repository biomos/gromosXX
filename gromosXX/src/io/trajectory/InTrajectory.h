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
   * it into simulation::System.
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
    template <math::boundary_enum b>
    InTrajectory & operator>>(simulation::System<b> &sys);
    /**
     * switch: require position (positionred) or not.
     */
    bool read_position;
    /**
     * switch: require velocity (velocityred) or not.
     */
    bool read_velocity;
    /**
     * switch: require box or not.
     */
    bool read_box;
    /**
     * switch: require box indices or not.
     */
    bool read_boxindices;
    
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
    template<math::boundary_enum b>
    bool _read_box(simulation::System<b> &sys, std::vector<std::string> &buffer);

    /**
     * read BOXINDICES block.
     */
    template<math::boundary_enum b>
    bool _read_boxindices(simulation::System<b> &sys,
			  std::vector<std::string> &buffer);
    
  };
  

} // io

// template and inline methods
#include "InTrajectory.tcc"

#endif
