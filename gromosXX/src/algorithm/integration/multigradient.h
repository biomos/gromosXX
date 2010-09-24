/**
 * @file multigradient.h
 * multiple gradients algorithm
 */

#ifndef INCLUDED_MULTIGRADIENT_H
#define INCLUDED_MULTIGRADIENT_H

namespace algorithm {

  /**
   * @class BezierPoint
   * helper class to represent a Bezier control point
   * @param time
   * @param value
   */
  class BezierPoint {
  public:

    BezierPoint(double time, double value) : time(time), value(value) {
    }

    BezierPoint(const BezierPoint & p) : time(p.time), value(p.value) {
    }

    BezierPoint & operator=(const BezierPoint & p) {
      time = p.time;
      value = p.value;
      return *this;
    }
    double time;
    double value;
  };

  /**
   * @class Bezier
   * A Bezier curve
   */
  class Bezier {
  public:
    std::string variable;
    double get_value(double time) const;
    std::vector<BezierPoint> control_points;
    std::string plot_ascii(double start_time, double end_time,
            unsigned int width, unsigned int height,
            std::string x_label = "", std::string y_label = "",
            const std::string & indent = "") const;
  };
  /**
   * @class Multi_Gradient
   * implements Multiple Gradients
   */
  class Multi_Gradient : public Algorithm
  {
    std::vector<Bezier> curves;
  public:
    /**
     * Constructor.
     */
    Multi_Gradient() : Algorithm("Multi_Gradient") {}

    /**
     * Destructor.
     */
    virtual ~Multi_Gradient() {}

    /**
     * apply the curves
     */
    virtual int apply(topology::Topology &topo,
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);
    /**
     * init
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

  };

} // algorithm

#endif
