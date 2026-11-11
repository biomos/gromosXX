/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file colvar_bias.cc
 * @brief Reusable bias/restraint model for collective variables.
 * This class converts an instantaneous collective-variable value into a
 * restrained value and evaluates the corresponding bias potential.
 *
 * Responsibilities:
 * - optional time averaging of the CV value
 * - optional transformed averaging, e.g. <r^-3>^-1/3 for distances
 * - full / upper-wall / lower-wall harmonic restraints
 * - optional linear-tail continuation
 * - optional periodic target wrapping, e.g. for dihedrals
 *
 * The Colvar itself only computes:
 *   s(x) and ds/dx
 *
 * This class computes:
 *   V(s), dV/ds
 *
 * The interaction layer then applies the chain rule:
 *   F = - dV/ds * ds/dx
 */


#include "../../stdheader.h"
#include "../../interaction/special/colvar/colvar_bias.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special



namespace interaction {
/**
 * @brief Default bias settings.
 *
 * By default this is an inactive/zero-strength two-sided harmonic restraint
 * without time averaging, periodicity, or linear-tail behavior.
 */
Colvar_Bias::Settings::Settings()
  : target(0.0),
    k(0.0),
    weight(1.0),
    linear_cutoff(0.0),
    rah(0),
    virial(0),
    wall_mode(wall_two_sided),
    periodic(false),
    period(0.0),
    paired_wall(false),
    wall_partner_target(0.0),
    explicit_periodic_branch(false),
    branch_delta(0.0),
    average_transform(average_none),
    force_scale_mode(forcescale_none),
    exponential_term(0.0),
    use_time_average(false),
    initialise_average_from_current(true)
{
}

Colvar_Bias::Result::Result()
  : value(0.0),
    energy(0.0),
    dE_dinstant(0.0),
    dE_dvalue(0.0),
    dE_dlinear_cutoff(0.0),
    dE_dk(0.0),
    force_scale(1.0),
    active(false),
    linear(false)
{
}

Colvar_Bias::Colvar_Bias()
  : m_settings(),
    m_average_initialised(false),
    m_average_state(0.0)
{
}

Colvar_Bias::Colvar_Bias(const Settings &settings)
  : m_settings(settings),
    m_average_initialised(false),
    m_average_state(0.0)
{
}

void Colvar_Bias::set_settings(const Settings &settings)
{
  m_settings = settings;
  m_average_initialised = false;
  m_average_state = 0.0;
}

void Colvar_Bias::update_settings_preserve_state(const Settings &settings)
{
  m_settings = settings;
}

const Colvar_Bias::Settings &Colvar_Bias::settings() const
{
  return m_settings;
}

/**
 * Averaging functionalities 
 */
bool Colvar_Bias::average_initialised() const
{
  return m_average_initialised;
}

double Colvar_Bias::average_state() const
{
  return m_average_state;
}

void Colvar_Bias::reset_average(double instant_value)
{
  m_average_state = transformed_value(instant_value);
  m_average_initialised = true;
}

bool Colvar_Bias::restore_average_state(double state, bool initialised)
{
  if (!std::isfinite(state)
      || (initialised && m_settings.average_transform == average_inverse_cubic && state <= 0.0))
    return false;
  m_average_state = state;
  m_average_initialised = initialised;
  return true;
}

double Colvar_Bias::transformed_value(double x) const
{
  if (m_settings.average_transform == average_inverse_cubic) {
    return std::pow(x, -3.0);
  }
  return x;
}

double Colvar_Bias::inverse_transformed_value(double y) const
{
  if (m_settings.average_transform == average_inverse_cubic) {
    return std::pow(y, -1.0 / 3.0);
  }
  return y;
}

double Colvar_Bias::d_transformed_dx(double x) const
{
  if (m_settings.average_transform == average_inverse_cubic) {
    return -3.0 * std::pow(x, -4.0);
  }
  return 1.0;
}

double Colvar_Bias::d_inverse_transformed_dy(double y) const
{
  if (m_settings.average_transform == average_inverse_cubic) {
    return -(1.0 / 3.0) * std::pow(y, -4.0 / 3.0);
  }
  return 1.0;
}

/**
 * @brief Difference between restrained value and target.
 *
 * For non-periodic CVs this is simply:
 *   value - target
 *
 * Full harmonic periodic biases use the shortest wrapped difference. Paired
 * walls instead share an image about the allowed arc's midpoint, so their
 * combined energy is continuous across the periodic branch.
 */
double Colvar_Bias::periodic_delta(double value, double target) const
{
  double delta = value - target;

  if (!m_settings.periodic || m_settings.period <= 0.0) {
    return delta;
  }

  if (m_settings.explicit_periodic_branch) {
    if (!std::isfinite(m_settings.branch_delta)) return std::numeric_limits<double>::quiet_NaN();
    // Match the legacy DELTA interval, including its endpoint convention.
    // Normalize DELTA while preserving the legacy distinction between zero
    // and a positive whole turn (the target lies at the opposite endpoint).
    double upper = std::fmod(m_settings.branch_delta, m_settings.period);
    if (upper < 0.0) upper += m_settings.period;
    if (upper == 0.0 && m_settings.branch_delta > 0.0) upper = m_settings.period;
    const double lower = upper - m_settings.period;
    if (delta < lower)
      delta += m_settings.period * std::ceil((lower - delta) / m_settings.period);
    if (delta > upper)
      delta -= m_settings.period * std::ceil((delta - upper) / m_settings.period);
    return delta;
  }

  if (m_settings.paired_wall && m_settings.wall_mode != wall_two_sided) {
    const bool lower = m_settings.wall_mode == wall_lower;
    const double lower_target = lower ? target : m_settings.wall_partner_target;
    const double upper_target = lower ? m_settings.wall_partner_target : target;
    double width = upper_target - lower_target;
    width -= m_settings.period * std::floor(width / m_settings.period);
    // Both walls use the same image about the allowed arc's midpoint. Their
    // sum is flat inside the arc and rises toward the opposite midpoint.
    delta = value - (lower_target + 0.5 * width);
    delta -= m_settings.period * std::floor(delta / m_settings.period + 0.5);
    return delta + (lower ? 0.5 * width : -0.5 * width);
  }

  delta -= m_settings.period * std::floor(delta / m_settings.period + 0.5);
  return delta;
}

/**
 * @brief Decide whether a one-sided restraint is active.
 *
 * wall_two_sided:
 *   active on both sides of the target
 *
 * wall_upper:
 *   active only if value > target
 *
 * wall_lower:
 *   active only if value < target
 */
bool Colvar_Bias::restraint_active(double delta) const
{
  if (m_settings.wall_mode == wall_upper && delta < 0.0) {
    return false;
  }
  if (m_settings.wall_mode == wall_lower && delta > 0.0) {
    return false;
  }
  return true;
}

/**
 * CORE functionality
 * @brief Evaluate bias energy and derivative for an instantaneous CV value.
 *
 * Main workflow:
 *
 * 1. Start from the instantaneous CV value s(t).
 * 2. Optionally update the time-averaged state.
 * 3. Convert the averaged state back to the restrained value.
 * 4. Compute delta = restrained_value - target.
 * 5. Apply wall logic.
 * 6. Evaluate harmonic or linear-tail potential.
 * 7. Return dE/ds_instant for force propagation.
 *
 * The interaction layer should multiply Result::dE_dinstant by the CV
 * derivative ds/dx.
 */
Colvar_Bias::Result Colvar_Bias::evaluate(double instant_value, bool advance_average)
{
  Result result;

  double restrained_value = instant_value;
  double dvalue_dinstant = 1.0;

  if (m_settings.use_time_average) {
    if (!m_average_initialised) {
      reset_average(instant_value);
    }

    if (advance_average) {
      const double current_transformed = transformed_value(instant_value);
      m_average_state = (1.0 - m_settings.exponential_term) * current_transformed
                      +        m_settings.exponential_term * m_average_state;
    }

    restrained_value = inverse_transformed_value(m_average_state);

    if (m_settings.force_scale_mode == forcescale_none) {
      dvalue_dinstant = 1.0;
    }
    else if (m_settings.force_scale_mode == forcescale_relaxation) {
      dvalue_dinstant = 1.0 - m_settings.exponential_term;
    }
    else {
      dvalue_dinstant = (1.0 - m_settings.exponential_term)
                      * d_inverse_transformed_dy(m_average_state)
                      * d_transformed_dx(instant_value);
    }
  }

  const double delta = periodic_delta(restrained_value, m_settings.target);

  result.value = restrained_value;
  result.force_scale = dvalue_dinstant;
  result.active = restraint_active(delta);

  if (!result.active) {
    return result;
  }

  const double abs_delta = std::fabs(delta);
  const bool use_linear_tail = (m_settings.linear_cutoff > 0.0
                             && abs_delta > m_settings.linear_cutoff);

  if (use_linear_tail) {
    result.linear = true;
    result.dE_dlinear_cutoff = m_settings.k * m_settings.weight
                            * (abs_delta - m_settings.linear_cutoff);
    result.dE_dk = m_settings.weight * m_settings.linear_cutoff
                 * (abs_delta - 0.5 * m_settings.linear_cutoff);

    result.dE_dvalue = m_settings.k * m_settings.weight
                     * m_settings.linear_cutoff
                     * (delta < 0.0 ? -1.0 : 1.0);
  }
  else {
    result.dE_dk = 0.5 * m_settings.weight * delta * delta;
    result.dE_dvalue = m_settings.k * m_settings.weight * delta;
  }

  result.energy = m_settings.k * result.dE_dk;
  result.dE_dinstant = result.dE_dvalue * dvalue_dinstant;
  return result;
}

Colvar_Bias::Lambda_Envelope Colvar_Bias::hidden_envelope(
    double lambda, unsigned int n, unsigned int m)
{
  // Combine the logarithms before exponentiation: separate factors can
  // produce 0*infinity even when their product is finite. Extended precision
  // also retains cancellation between the two derivative terms.
  const long double left = 2.0L * lambda;
  const long double right = 2.0L * (1.0L - lambda);
  const auto product = [=](unsigned int a, unsigned int b) -> long double {
    if ((a && left == 0) || (b && right == 0)) return 0.0L;
    const long double logarithm = (a ? a * std::log(std::fabs(left)) : 0.0L)
                                + (b ? b * std::log(std::fabs(right)) : 0.0L);
    // Integer powers remain meaningful just outside [0,1], too.
    const bool negative = (left < 0 && a % 2) != (right < 0 && b % 2);
    return (negative ? -1.0L : 1.0L) * std::exp(logarithm);
  };
  const long double derivative = (n ? 2.0L * n * product(n - 1, m) : 0.0L)
                               - (m ? 2.0L * m * product(n, m - 1) : 0.0L);
  return {static_cast<double>(product(n, m)), static_cast<double>(derivative)};
}

} // namespace interaction
