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
 * For periodic CVs, the shortest wrapped difference is used. This is needed
 * for angular CVs such as dihedrals.
 */
double Colvar_Bias::periodic_delta(double value, double target) const
{
  double delta = value - target;

  if (!m_settings.periodic || m_settings.period <= 0.0) {
    return delta;
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
Colvar_Bias::Result Colvar_Bias::evaluate(double instant_value)
{
  Result result;

  double restrained_value = instant_value;
  double dvalue_dinstant = 1.0;

  if (m_settings.use_time_average) {
    if (!m_average_initialised) {
      reset_average(instant_value);
    }

    const double old_average = m_average_state;
    const double current_transformed = transformed_value(instant_value);

    m_average_state = (1.0 - m_settings.exponential_term) * current_transformed
                    +        m_settings.exponential_term  * old_average;

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
    result.energy = m_settings.k * m_settings.weight * m_settings.linear_cutoff
                  * (abs_delta - 0.5 * m_settings.linear_cutoff);

    result.dE_dvalue = m_settings.k * m_settings.weight
                     * m_settings.linear_cutoff
                     * (delta < 0.0 ? -1.0 : 1.0);
  }
  else {
    result.energy = 0.5 * m_settings.k * m_settings.weight * delta * delta;
    result.dE_dvalue = m_settings.k * m_settings.weight * delta;
  }

  result.dE_dinstant = result.dE_dvalue * dvalue_dinstant;
  return result;
}

} // namespace interaction
