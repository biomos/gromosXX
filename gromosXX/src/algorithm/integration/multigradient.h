/**
 * @file multigradient.h
 * multiple gradients algorithm
 */

#ifndef INCLUDED_MULTIGRADIENT_H
#define INCLUDED_MULTIGRADIENT_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace algorithm {

  /**
   * @class Multi_Gradient
   * implements Multiple Gradients
   */
  class Multi_Gradient : public Algorithm {
  public:

    /**
     * @class ControlPoint
     * helper class to represent a control point
     * @param time
     * @param value
     */
    class ControlPoint {
    public:

      ControlPoint(double time, double value) : time(time), value(value) {
      }

      ControlPoint(const ControlPoint & p) : time(p.time), value(p.value) {
      }

      ControlPoint & operator=(const ControlPoint & p) {
        time = p.time;
        value = p.value;
        return *this;
      }
      double time;
      double value;
    };

    /**
     * @class Curve
     * an abstract base class for a curve
     */
    class Curve {
      std::vector<ControlPoint> m_control_points;
    public:
      virtual ~Curve() {}
      const std::vector<ControlPoint> & control_points() const {
        return m_control_points;
      }
      virtual void add_control_point(const ControlPoint & cp) {
        m_control_points.push_back(cp);
      }
      virtual void init() {}
      std::string variable;
      virtual double get_value(double time) const = 0;
      virtual void echo(std::ostream & os) const = 0;
      std::string plot_ascii(double start_time, double end_time,
              unsigned int width, unsigned int height,
              std::string x_label = "", std::string y_label = "",
              const std::string & indent = "") const;
    };

    /**
     * @class LinearInterpolation
     * points with linear interpolation between them
     */
    class LinearInterpolation : public Curve {
    public:
      double get_value(double time) const;
      void echo(std::ostream & os) const;
    };

     /**
     * @class SplineInterpolation
     * points with cubic spline interpolation between them
     */
    class SplineInterpolation : public Curve {
      // we need them as arrays for GSL...
      std::vector<double> times;
      std::vector<double> values;
      gsl_interp_accel *acc;
      gsl_spline *spline;
    public:
      SplineInterpolation() {
        acc = gsl_interp_accel_alloc();
        spline = NULL;
      }
      virtual ~SplineInterpolation() {
        if (spline != NULL)
          gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
      }
      void init();
      double get_value(double time) const;
      void echo(std::ostream & os) const;
    };

    /**
     * @class Bezier
     * A Bezier curve
     */
    class Bezier : public Curve {
    public:
      double get_value(double time) const;
      void echo(std::ostream & os) const;
    };

    /**
     * @class Oscillation
     * A cruve of type A sin(2pi/T*(t-DeltaT)) + b
     */
    class Oscillation : public Curve {
    public:
      double get_value(double time) const;
      void echo(std::ostream & os) const;
    };


    /**
     * Constructor.
     */
    Multi_Gradient() : Algorithm("Multi_Gradient") {
    }

    /**
     * Destructor.
     */
    virtual ~Multi_Gradient();

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

  private:
    std::vector<Curve*> curves;

  };

} // algorithm

#endif
