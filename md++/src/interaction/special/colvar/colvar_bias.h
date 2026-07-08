/**
 * @file colvar_bias.h
 * @brief Reusable bias/restraint model for collective variables.
 */

#ifndef INCLUDED_COLVAR_BIAS_H
#define INCLUDED_COLVAR_BIAS_H

#include <cmath>
#include <string>

namespace interaction {

  /**
   * Generic restraint/bias applied to a scalar collective variable.
   *
   * The bias owns state such as time averages.  The Colvar itself should stay
   * stateless except for its current value and derivatives.
   */
  class Colvar_Bias {
  public:
    enum Wall_Mode {
      wall_lower = -1,   // restrain only values below target
      wall_two_sided = 0,// restrain both sides of target
      wall_upper = 1     // restrain only values above target
    };

    enum Average_Transform {
      average_none = 0,
      average_identity = 1,
      average_inverse_cubic = 2
    };

    enum Force_Scale_Mode {
      forcescale_none = 0,
      forcescale_relaxation = 1,
      forcescale_chain_rule = 2
    };

    struct Settings {
      double target;
      double k;
      double weight;
      double linear_cutoff;
      int rah;
      unsigned int virial;
      Wall_Mode wall_mode;

      bool periodic;
      double period;

      Average_Transform average_transform;
      Force_Scale_Mode force_scale_mode;
      double exponential_term;
      bool use_time_average;
      bool initialise_average_from_current;

      Settings();
    };

    struct Result {
      double value;          // value actually restrained, possibly averaged
      double energy;
      double dE_dinstant;    // derivative wrt instantaneous colvar value
      double dE_dvalue;      // derivative wrt restrained value
      double force_scale;
      bool active;
      bool linear;

      Result();
    };

    Colvar_Bias();
    explicit Colvar_Bias(const Settings &settings);

    void set_settings(const Settings &settings);
    void update_settings_preserve_state(const Settings &settings);
    const Settings &settings() const;

    void reset_average(double instant_value);
    bool average_initialised() const;
    double average_state() const;

    Result evaluate(double instant_value);

  private:
    Settings m_settings;
    bool m_average_initialised;
    double m_average_state;

    double transformed_value(double x) const;
    double inverse_transformed_value(double y) const;
    double d_transformed_dx(double x) const;
    double d_inverse_transformed_dy(double y) const;

    double periodic_delta(double value, double target) const;
    bool restraint_active(double delta) const;
  };

} // namespace interaction

#endif
