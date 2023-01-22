/**
 * @file qm_link.cc
 * Implements methods of QM_Link
 */
#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../math/boundary_implementation.h"
#include "../../../math/periodicity.h"

#include "../../../util/debug.h"


#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

void interaction::QM_Link::update_cap_position(const math::Vec& qm_pos, const math::Vec& mm_pos, const double cap_length) const
  {
  math::Vec r_qmmm = mm_pos - qm_pos;
  DEBUG(15, "QM pos: " << math::v2s(qm_pos));
  DEBUG(15, "MM pos: " << math::v2s(mm_pos));
  DEBUG(15, "Link vector: " << math::v2s(r_qmmm));
  DEBUG(15, "Cap length: " << cap_length);
  
  this->pos = qm_pos + cap_length * (r_qmmm / math::abs(r_qmmm));
  DEBUG(15, "Capping atom position: " << math::v2s(this->pos));
}

void interaction::QM_Link::distribute_force(const math::Vec &qm_pos
                                          , const math::Vec &mm_pos
                                          , math::Vec &qm_force
                                          , math::Vec &mm_force) const
  {
  math::Vec r_mm_qm = mm_pos - qm_pos;
  
  // F_LQM = (1 - d_frac) * FL + d/d_QMMM^3 * (FL . r_MM-QM) * r_MM-QM
  // F_LQM = d_frac * FL - d/d_QMMM^3 * (FL . r_MM-QM) * r_MM-QM

  const double d_l_qm = math::abs(this->pos - qm_pos);
  const double d2_mm_qm = math::abs2(r_mm_qm);
  const double d_mm_qm = sqrt(d2_mm_qm);
  const double d3_mm_qm = d2_mm_qm * d_mm_qm;
  const double d_frac = d_l_qm / d_mm_qm;
  const double d3_frac = d_l_qm / d3_mm_qm;

  math::Vec l_qm_force = (1 - d_frac) * this->force;
  math::Vec l_mm_force = d_frac * this->force;

  const math::Vec b = d3_frac * math::dot(this->force, r_mm_qm) * r_mm_qm;

  l_qm_force += b;
  l_mm_force -= b;
  DEBUG(15, "Force on the capping atom: " << math::v2s(this->force));
  DEBUG(15, "Capping atom force distributed to QM atom: " << math::v2s(l_qm_force));
  DEBUG(15, "Capping atom force distributed to MM atom: " << math::v2s(l_mm_force));

  qm_force += l_qm_force;
  mm_force += l_mm_force;

  // Zero capping atom force, as it has been distributed
  this->force = 0.0;
}